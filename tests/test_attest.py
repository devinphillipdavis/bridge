"""bridge.attest: a manifest binds finished runs; a signature over it fails on any change.

The keys here are generated for the tests and sign only test fixtures."""
import json

import pytest
from osiris_governance import signing

from bridge import attest
from bridge import drift_check as dc
from dnalang.ledger import Ledger


def make_run(root, rows=3):
    """A minimal E2-shaped run: pre-registration, proof, ledger and telemetry."""
    (root / "PREREGISTRATION.md").write_text("plan\n")
    (root / "PREREGISTRATION.md.ots").write_bytes(b"proof")
    ledger = Ledger(root / "hardware" / "ledger.jsonl")
    ledger.append("preregistration", {"sha256": dc.sha256_file(root / "PREREGISTRATION.md")})
    for k in range(rows - 1):
        ledger.append("submit-intent", {"job": k})
    counts = [{"0": 990, "1": 10}, {"1": 985, "0": 15}, {"0": 900, "1": 100}]
    jobs = [{"window": (f"2026-10-26T15:{j:02d}:00.000000Z", f"2026-10-26T15:{j:02d}:03.500000Z"),
             "counts": counts} for j in range(3)]
    lines = dc.write_telemetry(jobs, [3], 1000, "attest-test", {"kind": "HARDWARE", "node": "fake",
                                                                "generator": "test"}, "MEASURED")
    (root / "hardware" / "telemetry.jsonl").write_bytes(b"".join(line + b"\n" for line in lines))
    (root / "hardware" / "analysis.json").write_text("{}\n")
    return root / "ATTESTATION.json"


@pytest.fixture
def keys():
    private = signing.generate_private_key_pem()
    return private, signing.public_key_pem(private)


def test_build_binds_the_ledger_head_telemetry_and_every_file(tmp_path):
    manifest = attest.build(make_run(tmp_path, rows=4), ["E2"])
    run = manifest["runs"]["E2"]
    rows = list(Ledger(tmp_path / "hardware" / "ledger.jsonl"))
    assert run["ledger"]["rows"] == 4 and run["ledger"]["head"] == rows[-1]["hash"]
    assert run["telemetry"]["records"] == 3 and run["telemetry"]["run_id"] == "attest-test"
    assert set(run["files"]) == {"hardware/ledger.jsonl", "hardware/telemetry.jsonl", "hardware/analysis.json"}
    assert run["preregistration"] == {"path": "PREREGISTRATION.md", "ots": "PREREGISTRATION.md.ots",
                                      "sha256": dc.sha256_file(tmp_path / "PREREGISTRATION.md")}
    assert manifest["does_not_show"] and "Retrospective" in manifest["statement"]
    assert attest.check(tmp_path / "ATTESTATION.json") == []


def test_build_refuses_a_broken_chain_and_never_rewrites(tmp_path):
    out = make_run(tmp_path)
    path = tmp_path / "hardware" / "ledger.jsonl"
    lines = path.read_text().splitlines()
    row = json.loads(lines[1])
    row["job"] = 99
    path.write_text("\n".join([lines[0], json.dumps(row)] + lines[2:]) + "\n")
    with pytest.raises(ValueError, match="does not verify"):
        attest.build(out, ["E2"])
    out.write_text("{}")
    with pytest.raises(FileExistsError):
        attest.build(out, ["E2"])


def test_build_refuses_altered_telemetry_and_an_empty_ledger(tmp_path):
    out = make_run(tmp_path)
    path = tmp_path / "hardware" / "telemetry.jsonl"
    path.write_bytes(path.read_bytes().replace(b'"attest-test"', b'"attest-other"', 1))
    with pytest.raises(ValueError, match="telemetry does not verify"):
        attest.build(out, ["E2"])
    (tmp_path / "hardware" / "ledger.jsonl").write_text("")
    with pytest.raises(ValueError, match="ledger is empty"):
        attest.build(out, ["E2"])


def test_signed_manifest_verifies(tmp_path, keys):
    out = make_run(tmp_path)
    attest.build(out, ["E2"])
    record = attest.sign(out, keys[0])
    assert record["manifest_sha256"] == dc.sha256_file(out)
    assert "PRIVATE" not in attest.signature_path(out).read_text()
    assert attest.verify(out, keys[1]) == []


def test_verify_fails_on_a_changed_or_extra_file(tmp_path, keys):
    out = make_run(tmp_path)
    attest.build(out, ["E2"])
    attest.sign(out, keys[0])
    (tmp_path / "hardware" / "extra.json").write_text("{}")
    assert attest.verify(out, keys[1]) == ["E2: extra hardware/extra.json"]
    (tmp_path / "hardware" / "extra.json").unlink()
    (tmp_path / "hardware" / "analysis.json").unlink()
    assert attest.verify(out, keys[1]) == ["E2: missing hardware/analysis.json"]
    Ledger(tmp_path / "hardware" / "ledger.jsonl").append("result", {"job": 0})
    problems = attest.verify(out, keys[1])
    assert "E2: ledger differs" in problems and "E2: changed hardware/ledger.jsonl" in problems


def test_verify_fails_on_a_truncated_manifest_or_another_key(tmp_path, keys):
    out = make_run(tmp_path)
    attest.build(out, ["E2"])
    attest.sign(out, keys[0])
    other = signing.public_key_pem(signing.generate_private_key_pem())
    assert attest.verify(out, other) == ["the signature names a different key",
                                         "the signature does not match the manifest"]
    out.write_bytes(out.read_bytes()[:-2] + b"\n")
    assert "the signature does not match the manifest" in attest.verify(out, keys[1])


def test_sign_refuses_changed_files_and_a_second_signature(tmp_path, keys):
    out = make_run(tmp_path)
    attest.build(out, ["E2"])
    (tmp_path / "PREREGISTRATION.md").write_text("another plan\n")
    with pytest.raises(ValueError, match="preregistration differs"):
        attest.sign(out, keys[0])
    (tmp_path / "PREREGISTRATION.md").write_text("plan\n")
    attest.sign(out, keys[0])
    with pytest.raises(FileExistsError):
        attest.sign(out, keys[0])
