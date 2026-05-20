"""Shared visualization helpers."""


def has_cds(exon: dict) -> bool:
    """Return True when exon has a valid CDS interval."""
    return (
        exon.get("cds_start") is not None
        and exon.get("cds_end") is not None
        and exon["cds_end"] > exon["cds_start"]
    )
