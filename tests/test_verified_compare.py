"""bridge.verified_compare: decisions are recomputed from raw counts (needs [governance])."""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
pytest.importorskip("qiskit_aer")
pytest.importorskip("osiris_governance")

from osiris_governance.contracts import PermitAuthority  # noqa: E402
from osiris_governance.ledger import DynamicEvidenceLedger  # noqa: E402
from osiris_governance.models.ledger import LedgerEventType  # noqa: E402

import bridge.verified_compare as vc  # noqa: E402
from bridge.compare import HARD, compare  # noqa: E402

AUTHORITY = PermitAuthority(b"k" * 32)
SMALL = dict(seeds=[0], budget=10, cal=HARD, verify_shots=128)


def test_seals_raw_evidence_and_blocks_single_seed(tmp_path):
    ledger = DynamicEvidenceLedger(ledger_id="test-bridge-ledger", permit_authority=AUTHORITY)
    out = vc.verified_compare(workdir=tmp_path, ledger=ledger, permit_authority=AUTHORITY, **SMALL)
    assert out["verified_evidence_sealed"] is True
    assert all(i["matched"] for i in out["evidence_integrity"])
    for i, row in zip(out["evidence_integrity"], [r[m] for r in out["rows"] for m in ("organism", "ga")]):
        assert i["recomputed"] == pytest.approx(row["verify"], abs=1e-12)
    decision = out["statistical_promotion_decision"]
    assert decision["promoted"] is False
    assert decision["reason_code"] == "PSEUDO_REPLICATION_BLOCKED"
    assert decision["search_evaluations_checked"] > 0
    verified = [e for e in ledger.events if e.event_type == LedgerEventType.RAW_DATA_VERIFIED]
    assert len(verified) == 2
    assert (tmp_path / "verified_compare.json").exists()


def test_tampered_score_caught_from_counts(tmp_path, monkeypatch):
    def inflated(**kw):
        res = compare(**kw)
        res["rows"][0]["organism"]["verify"] += 0.01
        return res

    monkeypatch.setattr(vc, "compare", inflated)
    out = vc.verified_compare(workdir=tmp_path, **SMALL)
    assert out["verified_evidence_sealed"] is False
    assert out["statistical_promotion_decision"]["promoted"] is False
    assert out["statistical_promotion_decision"]["reason_code"] == "EVIDENCE_INTEGRITY_MISMATCH"


def test_ledger_requires_permit_authority(tmp_path):
    with pytest.raises(ValueError, match="permit_authority"):
        vc.verified_compare(workdir=tmp_path, ledger=DynamicEvidenceLedger(), **SMALL)
