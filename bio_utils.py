"""
bio_utils.py — BioPython helpers for parsing and summarising sequence uploads.

Keeps all biology-specific logic separate from the Streamlit UI layer.
"""

from __future__ import annotations

import io
import statistics
from dataclasses import dataclass, field
from typing import List, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@dataclass
class SequenceSummary:
    """Container for parsed FASTA summary statistics."""

    file_name: str
    format: str  # e.g. "Protein FASTA (.faa)"
    sequence_count: int
    total_residues: int
    avg_length: float
    min_length: int
    max_length: int
    sample_ids: List[str] = field(default_factory=list)  # first few IDs

    @property
    def size_label(self) -> str:
        """Human-friendly total residue count."""
        if self.total_residues >= 1_000_000:
            return f"{self.total_residues / 1_000_000:.1f}M residues"
        if self.total_residues >= 1_000:
            return f"{self.total_residues / 1_000:.1f}K residues"
        return f"{self.total_residues} residues"


def parse_fasta(file_obj: io.IOBase, file_name: str) -> SequenceSummary:
    """
    Parse an uploaded FASTA file and return a compact summary.

    Parameters
    ----------
    file_obj : file-like
        An open binary or text stream (e.g. Streamlit UploadedFile).
    file_name : str
        Original filename — used to detect format and stored in summary.

    Returns
    -------
    SequenceSummary

    Raises
    ------
    ValueError
        If the file contains no valid FASTA sequences.
    """
    # Streamlit UploadedFile is bytes-mode; SeqIO needs a text handle.
    text_handle = io.TextIOWrapper(file_obj, encoding="utf-8")

    records: List[SeqRecord] = list(SeqIO.parse(text_handle, "fasta"))

    if not records:
        raise ValueError(
            "No valid FASTA sequences found. "
            "Please upload a file in FASTA format (.faa or .fasta)."
        )

    lengths = [len(r.seq) for r in records]

    # Detect if protein vs nucleotide (simple heuristic)
    # Protein alphabets commonly include letters outside {A,C,G,T,N,U}
    sample_seq = str(records[0].seq).upper()[:200]
    nuc_chars = set("ACGTNU")
    non_nuc_ratio = sum(1 for c in sample_seq if c not in nuc_chars) / max(len(sample_seq), 1)
    is_protein = non_nuc_ratio > 0.1

    ext = file_name.rsplit(".", 1)[-1].lower() if "." in file_name else ""
    if ext in ("faa",):
        fmt_label = "Protein FASTA (.faa)"
    elif ext in ("fna", "ffn"):
        fmt_label = "Nucleotide FASTA (.{})".format(ext)
    elif is_protein:
        fmt_label = "Protein FASTA"
    else:
        fmt_label = "Nucleotide FASTA"

    return SequenceSummary(
        file_name=file_name,
        format=fmt_label,
        sequence_count=len(records),
        total_residues=sum(lengths),
        avg_length=round(statistics.mean(lengths), 1),
        min_length=min(lengths),
        max_length=max(lengths),
        sample_ids=[r.id for r in records[:5]],
    )


def generate_example_fasta() -> tuple[bytes, str]:
    """
    Return a small synthetic protein FASTA as (bytes, filename).

    Useful for the "Load Example Dataset" button so the parsing
    pipeline runs on real data even in demo mode.
    """
    file_name = "E_coli_K12_proteome_sample.faa"

    # 12 realistic-looking E. coli protein entries (truncated sequences)
    entries = [
        ("thrA", "Bifunctional aspartokinase/homoserine dehydrogenase I",
         "MRVLKFGGTSVANAERFLRVAKDLAENNRGARVLVVCSEITAVTFRGPSETH"
         "LDSMVGQALFGDGAGAVIVGSDPDLSVERPLYAEITLLATQEASKQRSF"),
        ("thrB", "Homoserine kinase",
         "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNLKQ"
         "ERFLNRVTQVSIEQPKSTEELVKRILQAGEFTWVIMDPHTDAKAIAELSDK"),
        ("thrC", "Threonine synthase",
         "MKLYNLKDHNEQVSFAQAVTQGLGKNQGLFFPHDLPEFSLTEIDEMLKLD"
         "FAKLAGAKNFLPKERTFDQPELEALAQGARALCADAVLVDSATAALKRGI"),
        ("dnaK", "Chaperone protein DnaK",
         "MGKIIGIDLGTTNSCVAIMDGTTPRVLENAEGDRTTPSIIAYTQDGETLVG"
         "QPAKRQAVTNPQNTLFAIKRLIGRRFQDEEVQRDVSIMPFKIIAADNGDAW"),
        ("groEL", "60 kDa chaperonin",
         "MAAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITK"
         "DGVSVAREIELEDKFENMGAQMVKEVASKANDAAGDGTTTATVLAQAIIAEL"),
        ("rpoB", "DNA-directed RNA polymerase subunit beta",
         "MVYSYTEKKRIRKDFGKRPQVLDIPYLLSIFLSKVHIPAKKQAPALLGKNIK"
         "DVITLREFEDRIIAQFLKTNSHRSGALAALCALVARECGTVPNUHATKHEI"),
        ("gapA", "Glyceraldehyde-3-phosphate dehydrogenase A",
         "MTIKVGINGFGRIGRIVFRAAQKRSDIEIVAINDLLDADYMAYMLKYDSTHG"
         "RFDGTVEVKDGHLIVNGKKIRVTAERDPANLKWDEVGVDVVAEATGQFHSE"),
        ("ompA", "Outer membrane protein A",
         "MKKTAIAIAVALAGFATVAAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGP"
         "THENQLGAGAFGGYQVNPYVGFEMGYDWLGRMPYKGSVENGAYKAQGVQLT"),
        ("atpA", "ATP synthase subunit alpha",
         "MQLNSTEISALQLKVKELEQHEHLKKAIAQIDKFREVGAADTSKLTDAMADI"
         "LGKIDHQLLERLEAAQQSDSQILAIGLKKLASMTPFVSQNSSYLAERYGLVV"),
        ("ftsZ", "Cell division protein FtsZ",
         "MFEPMELTNDAVIKVIGVGGGGGNAVEHMVRERIEGPADLVNVDLAAQALRK"
         "EKDPSVLSQMGADLNLVFLTGQTGAGKSTTSSAAITKVLREIAKEREVPIL"),
        ("recA", "Protein RecA",
         "MAIDENKQKALAAALGQIEKQFGKGSIMRLGEDRSMDVETISTGSLSLDIALGAG"
         "GLPMGRIVEIYGPESSGKTTLTLQVIAAAQREGKTCAFIDAEHALDPIYA"),
        ("tufA", "Elongation factor Tu 1",
         "MSKEKFERTKPHVNVGTIGHVDHGKTTLTAAITTVLAKTYGGAARAFDQIDNA"
         "PEEKARGITINTSHVEYDTPTRHYAHVDCPGHADYVKNMITGAAQMDGAIL"),
    ]

    lines: list[str] = []
    for gene, desc, seq in entries:
        lines.append(f">sp|P0A{gene[:3].upper()}0|{gene.upper()}_ECOLI {desc} OS=Escherichia coli (strain K12)")
        # Wrap sequence at 70 chars
        for i in range(0, len(seq), 70):
            lines.append(seq[i : i + 70])

    content = "\n".join(lines) + "\n"
    return content.encode("utf-8"), file_name
