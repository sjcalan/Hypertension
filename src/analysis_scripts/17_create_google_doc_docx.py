#!/usr/bin/env python3
"""Create a Google-Docs-ready DOCX from the current LaTeX manuscript.

This avoids external converters so it can run on the server with only the
Python standard library. The output is a .docx file that Google Docs can import.
"""

from __future__ import annotations

import html
import re
import struct
import zipfile
from dataclasses import dataclass
from pathlib import Path


BASE_DIR = Path("/userhome/cs3/u3011656/hypertension/hypertension")
PROJECT_DIR = BASE_DIR / "analysis" / "overleaf_hypertension_visit_trajectory_2026-05-25"
OUT_DIR = BASE_DIR / "analysis" / "outputs" / "google_doc_hypertension_visit_trajectory_2026-06-11"
OUT_DOCX = OUT_DIR / "hypertension_visit_trajectory_google_doc_2026-06-11.docx"

LABEL_MAP = {
    "fig:flow": "Figure 1",
    "fig:retention": "Figure 2",
    "fig:clinical_trajectories": "Figure 3",
    "fig:journey_appendix": "Supplementary Figure S1",
    "fig:heatmap_appendix": "Supplementary Figure S2",
    "tab:baseline": "Table 1",
    "tab:cluster_comparison": "Table 2",
    "tab:model_selection": "Supplementary Table S1",
    "tab:excluded_modeled": "Supplementary Table S2",
    "tab:administrative_followup": "Supplementary Table S3",
    "tab:bp_opportunity": "Supplementary Table S4",
    "tab:gbtm_diagnostics_status": "Supplementary Table S5",
}


@dataclass
class BibEntry:
    key: str
    authors_text: str
    cite_text: str
    year: str
    title: str
    journal: str
    volume: str = ""
    number: str = ""
    pages: str = ""
    doi: str = ""

    def reference_text(self) -> str:
        bits = [f"{self.authors_text} ({self.year}). {self.title}. {self.journal}"]
        tail = ""
        if self.volume:
            tail += self.volume
        if self.number:
            tail += f"({self.number})"
        if self.pages:
            tail += f": {self.pages}"
        if tail:
            bits[-1] += f", {tail}"
        if self.doi:
            bits[-1] += f". https://doi.org/{self.doi}"
        else:
            bits[-1] += "."
        return bits[-1]


def x(text: str) -> str:
    return html.escape(text, quote=True)


def parse_bib(path: Path) -> dict[str, BibEntry]:
    raw = path.read_text()
    entries: dict[str, BibEntry] = {}
    for match in re.finditer(r"@article\{([^,]+),(.*?)\n\}", raw, flags=re.S):
        key, body = match.group(1).strip(), match.group(2)
        fields = {}
        for fm in re.finditer(r"(\w+)\s*=\s*\{(.*?)\}\s*,?", body, flags=re.S):
            fields[fm.group(1).lower()] = re.sub(r"\s+", " ", fm.group(2).strip())
        authors = fields.get("author", "")
        author_parts = [a.strip() for a in authors.split(" and ")]
        first_author_last = author_parts[0].split(",")[0] if author_parts else key
        cite = f"{first_author_last} et al., {fields.get('year', '')}"
        display_authors = []
        for author in author_parts:
            if author.lower() == "others":
                display_authors.append("et al.")
            elif "," in author:
                last, rest = [p.strip() for p in author.split(",", 1)]
                initials = " ".join(part[:1] for part in re.split(r"\s+", rest) if part)
                display_authors.append(f"{last} {initials}".strip())
            else:
                display_authors.append(author)
        entries[key] = BibEntry(
            key=key,
            authors_text=", ".join(display_authors),
            cite_text=cite,
            year=fields.get("year", ""),
            title=fields.get("title", ""),
            journal=fields.get("journal", ""),
            volume=fields.get("volume", ""),
            number=fields.get("number", ""),
            pages=fields.get("pages", "").replace("--", "-"),
            doi=fields.get("doi", ""),
        )
    return entries


BIB = parse_bib(PROJECT_DIR / "references.bib")


def citation_repl(match: re.Match[str]) -> str:
    keys = [k.strip() for k in match.group(1).split(",")]
    cites = [BIB[k].cite_text if k in BIB else k for k in keys]
    return "(" + "; ".join(cites) + ")"


def strip_latex_comments(text: str) -> str:
    """Remove unescaped LaTeX comments without treating escaped \\% as comments."""
    cleaned_lines = []
    for line in text.splitlines():
        search_from = 0
        while True:
            percent_idx = line.find("%", search_from)
            if percent_idx == -1:
                cleaned_lines.append(line)
                break

            backslash_count = 0
            cursor = percent_idx - 1
            while cursor >= 0 and line[cursor] == "\\":
                backslash_count += 1
                cursor -= 1

            if backslash_count % 2 == 0:
                cleaned_lines.append(line[:percent_idx])
                break

            search_from = percent_idx + 1
    return "\n".join(cleaned_lines)


def latex_to_text(text: str) -> str:
    text = strip_latex_comments(text.strip())
    text = re.sub(r"\\citep\{([^}]+)\}", citation_repl, text)
    for label, replacement in LABEL_MAP.items():
        text = text.replace(rf"\ref{{{label}}}", replacement)
    text = text.replace("~", " ")
    text = text.replace(r"\bp{}", "blood pressure")
    text = text.replace(r"\bp", "blood pressure")
    text = text.replace(r"\%", "%")
    text = text.replace(r"\&", "&")
    text = text.replace(r"\#", "#")
    text = text.replace(r"\$", "$")
    text = text.replace(r"\_", "_")
    text = text.replace(r"\pm", "±")
    text = text.replace(r"\newline", " ")
    text = text.replace(r"``", '"').replace(r"''", '"')
    text = re.sub(r"\\texttt\{([^{}]*)\}", r"\1", text)
    text = re.sub(r"\\textit\{([^{}]*)\}", r"\1", text)
    text = re.sub(r"\\emph\{([^{}]*)\}", r"\1", text)
    text = re.sub(r"\\textbf\{([^{}]*)\}", r"\1", text)
    text = re.sub(r"\$(.*?)\$", lambda m: m.group(1), text)
    text = re.sub(r"\\[a-zA-Z]+\*?(?:\[[^\]]*\])?(?:\{([^{}]*)\})?", lambda m: m.group(1) or "", text)
    text = text.replace("{", "").replace("}", "")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def png_size(path: Path) -> tuple[int, int]:
    with path.open("rb") as f:
        header = f.read(24)
    if header[:8] != b"\x89PNG\r\n\x1a\n":
        return (1200, 800)
    width, height = struct.unpack(">II", header[16:24])
    return width, height


class DocxBuilder:
    def __init__(self) -> None:
        self.body: list[str] = []
        self.rels: list[tuple[str, str, str]] = []
        self.media: list[tuple[Path, str]] = []
        self.next_rid = 1
        self.next_pic_id = 1

    def paragraph(self, text: str = "", style: str | None = None, *, bold: bool = False, italic: bool = False, align: str | None = None) -> None:
        text = text or ""
        ppr = ""
        if style:
            ppr += f'<w:pStyle w:val="{style}"/>'
        if align:
            ppr += f'<w:jc w:val="{align}"/>'
        ppr_xml = f"<w:pPr>{ppr}</w:pPr>" if ppr else ""
        rpr = ""
        if bold:
            rpr += "<w:b/>"
        if italic:
            rpr += "<w:i/>"
        rpr_xml = f"<w:rPr>{rpr}</w:rPr>" if rpr else ""
        runs = []
        parts = text.split("\n")
        for i, part in enumerate(parts):
            if i:
                runs.append("<w:r><w:br/></w:r>")
            runs.append(f'<w:r>{rpr_xml}<w:t xml:space="preserve">{x(part)}</w:t></w:r>')
        self.body.append(f"<w:p>{ppr_xml}{''.join(runs)}</w:p>")

    def heading(self, text: str, level: int) -> None:
        style = "Heading1" if level == 1 else "Heading2"
        self.paragraph(text, style=style)

    def add_table(self, rows: list[list[str]]) -> None:
        if not rows:
            return
        ncols = max(len(row) for row in rows)
        tbl_rows = []
        for row_idx, row in enumerate(rows):
            cells = []
            padded = row + [""] * (ncols - len(row))
            for cell in padded:
                shade = '<w:shd w:fill="EDEDED"/>' if row_idx == 0 else ""
                tcpr = f'<w:tcPr><w:tcW w:w="{int(9000 / max(ncols, 1))}" w:type="dxa"/>{shade}</w:tcPr>'
                cells.append(
                    f"<w:tc>{tcpr}<w:p><w:r><w:t>{x(cell)}</w:t></w:r></w:p></w:tc>"
                )
            tbl_rows.append("<w:tr>" + "".join(cells) + "</w:tr>")
        borders = (
            '<w:tblBorders><w:top w:val="single" w:sz="4" w:space="0" w:color="999999"/>'
            '<w:left w:val="single" w:sz="4" w:space="0" w:color="999999"/>'
            '<w:bottom w:val="single" w:sz="4" w:space="0" w:color="999999"/>'
            '<w:right w:val="single" w:sz="4" w:space="0" w:color="999999"/>'
            '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="CCCCCC"/>'
            '<w:insideV w:val="single" w:sz="4" w:space="0" w:color="CCCCCC"/></w:tblBorders>'
        )
        self.body.append(
            '<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/>'
            + borders
            + "</w:tblPr>"
            + "".join(tbl_rows)
            + "</w:tbl>"
        )
        self.paragraph("")

    def add_image(self, path: Path, caption: str | None = None, max_width_in: float = 6.2) -> None:
        if not path.exists():
            self.paragraph(f"[Missing figure: {path.name}]", italic=True)
            return
        rid = f"rId{self.next_rid}"
        self.next_rid += 1
        media_name = f"image{len(self.media) + 1}{path.suffix.lower()}"
        self.media.append((path, f"word/media/{media_name}"))
        self.rels.append((rid, "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image", f"media/{media_name}"))

        width_px, height_px = png_size(path) if path.suffix.lower() == ".png" else (1200, 800)
        width_in = min(max_width_in, width_px / 300.0)
        if width_in < 4.0:
            width_in = min(max_width_in, width_px / 150.0)
        height_in = width_in * height_px / width_px
        cx, cy = int(width_in * 914400), int(height_in * 914400)
        pic_id = self.next_pic_id
        self.next_pic_id += 1
        drawing = f"""
<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:drawing>
<wp:inline distT="0" distB="0" distL="0" distR="0">
<wp:extent cx="{cx}" cy="{cy}"/>
<wp:docPr id="{pic_id}" name="Picture {pic_id}"/>
<wp:cNvGraphicFramePr><a:graphicFrameLocks noChangeAspect="1"/></wp:cNvGraphicFramePr>
<a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">
<pic:pic><pic:nvPicPr><pic:cNvPr id="{pic_id}" name="{x(path.name)}"/><pic:cNvPicPr/></pic:nvPicPr>
<pic:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>
<pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{cx}" cy="{cy}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>
</pic:pic></a:graphicData></a:graphic>
</wp:inline></w:drawing></w:r></w:p>
"""
        self.body.append(drawing)
        if caption:
            self.paragraph(caption, style="Caption", italic=True)

    def save(self, path: Path) -> None:
        section = '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="720" w:right="720" w:bottom="720" w:left="720" w:header="720" w:footer="720" w:gutter="0"/></w:sectPr>'
        document_xml = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:wpc="http://schemas.microsoft.com/office/word/2010/wordprocessingCanvas"
 xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006"
 xmlns:o="urn:schemas-microsoft-com:office:office"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
 xmlns:m="http://schemas.openxmlformats.org/officeDocument/2006/math"
 xmlns:v="urn:schemas-microsoft-com:vml"
 xmlns:wp14="http://schemas.microsoft.com/office/word/2010/wordprocessingDrawing"
 xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
 xmlns:w10="urn:schemas-microsoft-com:office:word"
 xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"
 xmlns:w14="http://schemas.microsoft.com/office/word/2010/wordml"
 xmlns:wpg="http://schemas.microsoft.com/office/word/2010/wordprocessingGroup"
 xmlns:wpi="http://schemas.microsoft.com/office/word/2010/wordprocessingInk"
 xmlns:wne="http://schemas.microsoft.com/office/word/2006/wordml"
 xmlns:wps="http://schemas.microsoft.com/office/word/2010/wordprocessingShape"
 xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"
 xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture"
 mc:Ignorable="w14 wp14"><w:body>{''.join(self.body)}{section}</w:body></w:document>'''
        rels_xml = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        for rid, rel_type, target in self.rels:
            rels_xml += f'<Relationship Id="{rid}" Type="{rel_type}" Target="{target}"/>'
        rels_xml += "</Relationships>"
        styles_xml = styles()
        content_types = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Default Extension="png" ContentType="image/png"/>
<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>'''
        package_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>'''
        core = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
 xmlns:dc="http://purl.org/dc/elements/1.1/"
 xmlns:dcterms="http://purl.org/dc/terms/"
 xmlns:dcmitype="http://purl.org/dc/dcmitype/"
 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<dc:title>Hypertension Care Retention Trajectories and Blood Pressure Outcomes</dc:title>
<dc:creator>Jiacheng Song; Jiancheng Ye</dc:creator>
</cp:coreProperties>'''
        app = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"
 xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
<Application>Codex DOCX Builder</Application></Properties>'''
        path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("[Content_Types].xml", content_types)
            z.writestr("_rels/.rels", package_rels)
            z.writestr("word/document.xml", document_xml)
            z.writestr("word/_rels/document.xml.rels", rels_xml)
            z.writestr("word/styles.xml", styles_xml)
            z.writestr("docProps/core.xml", core)
            z.writestr("docProps/app.xml", app)
            for src, dst in self.media:
                z.write(src, dst)


def styles() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:qFormat/><w:rPr><w:sz w:val="22"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:qFormat/><w:pPr><w:jc w:val="center"/></w:pPr><w:rPr><w:b/><w:sz w:val="32"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Subtitle"><w:name w:val="Subtitle"/><w:qFormat/><w:pPr><w:jc w:val="center"/></w:pPr><w:rPr><w:sz w:val="22"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:qFormat/><w:pPr><w:spacing w:before="240" w:after="120"/></w:pPr><w:rPr><w:b/><w:sz w:val="28"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/><w:qFormat/><w:pPr><w:spacing w:before="180" w:after="80"/></w:pPr><w:rPr><w:b/><w:sz w:val="24"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Caption"><w:name w:val="Caption"/><w:qFormat/><w:rPr><w:i/><w:sz w:val="20"/></w:rPr></w:style>
</w:styles>'''


def table_rows_from_tex(path: Path) -> tuple[str, list[list[str]], str]:
    text = path.read_text()
    caption = latex_to_text(re.search(r"\\caption\{(.+?)\}", text, flags=re.S).group(1))
    note_match = re.search(r"\\item\s+(.+?)\n\\end\{tablenotes\}", text, flags=re.S)
    note = latex_to_text(note_match.group(1)) if note_match else ""
    body_match = re.search(r"\\toprule\n(.+?)\n\\bottomrule", text, flags=re.S)
    rows: list[list[str]] = []
    if body_match:
        for raw_line in body_match.group(1).splitlines():
            line = raw_line.strip()
            if not line or line in {r"\midrule"}:
                continue
            line = re.sub(r"\\\\\[.*?\]", r"\\", line)
            section = re.search(r"\\multicolumn\{\d+\}\{[^}]+\}\{\\textit\{(.+?)\}\}", line)
            if section:
                rows.append([latex_to_text(section.group(1))])
                continue
            if line.startswith(r"\addlinespace"):
                section = re.search(r"\\textit\{(.+?)\}", line)
                if section:
                    rows.append([latex_to_text(section.group(1))])
                continue
            line = re.sub(r"\\\\\s*$", "", line)
            cells = [latex_to_text(c) for c in line.split(" & ")]
            if len(cells) > 1:
                rows.append(cells)
    return caption, rows, note


def figure_caption(block: str) -> str:
    matches = re.findall(r"\\caption\{(.+?)\}", block, flags=re.S)
    return latex_to_text(matches[-1]) if matches else ""


def add_flow_figure(doc: DocxBuilder, caption: str) -> None:
    doc.paragraph("Figure 1. Cohort flow", style="Caption", italic=True)
    doc.add_table(
        [
            ["Step", "Patients"],
            ["Visit-level hypertension EHR cohort", "26,541"],
            ["Trajectory-eligible cohort: at least 3 visits", "17,857"],
            ["Relative-time visit panel: modeled cohort", "17,802"],
            ["Cluster 1: Long-span sparse", "2,703"],
            ["Cluster 2: Sustained moderate", "5,081"],
            ["Cluster 3: Short-span intensive early", "3,404"],
            ["Cluster 4: Frequent intermediate", "6,614"],
        ]
    )
    if caption:
        doc.paragraph(caption, style="Caption", italic=True)


def process_figure_block(doc: DocxBuilder, block: str) -> None:
    caption = figure_caption(block)
    if "tikzpicture" in block:
        add_flow_figure(doc, caption)
        return
    include_paths = re.findall(r"\\includegraphics(?:\[[^\]]+\])?\{([^}]+)\}", block)
    subcaptions = re.findall(r"\\caption\{(.+?)\}", block, flags=re.S)
    outer_caption = latex_to_text(subcaptions[-1]) if subcaptions else caption
    single_subcaps = [latex_to_text(c) for c in subcaptions[:-1]]
    for i, inc in enumerate(include_paths):
        img = PROJECT_DIR / "figures" / inc
        if not img.exists() and not Path(inc).suffix:
            img = PROJECT_DIR / "figures" / f"{inc}.png"
        cap = single_subcaps[i] if i < len(single_subcaps) else None
        doc.add_image(img, cap, max_width_in=5.9)
    if outer_caption:
        doc.paragraph(outer_caption, style="Caption", italic=True)


def parse_main() -> None:
    text = (PROJECT_DIR / "main.tex").read_text()
    title = latex_to_text(re.search(r"\\title\{(.+?)\}", text, flags=re.S).group(1))
    author_raw = re.search(r"\\author\{(.+?)\}", text, flags=re.S).group(1).replace(r"\and", ";")
    author = latex_to_text(author_raw)
    author = re.sub(r"\s*;\s*", "; ", author)
    abstract = latex_to_text(re.search(r"\\begin\{abstract\}(.+?)\\end\{abstract\}", text, flags=re.S).group(1))
    content = text.split(r"\end{abstract}", 1)[1].split(r"\end{document}", 1)[0]

    doc = DocxBuilder()
    doc.paragraph(title, style="Title")
    doc.paragraph(author, style="Subtitle")
    doc.heading("Abstract", 1)
    doc.paragraph(abstract)

    paragraph_buf: list[str] = []
    refs_added = False

    def flush_paragraph() -> None:
        nonlocal paragraph_buf
        if paragraph_buf:
            txt = latex_to_text(" ".join(paragraph_buf))
            if txt:
                doc.paragraph(txt)
            paragraph_buf = []

    lines = content.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            flush_paragraph()
            i += 1
            continue
        if line.startswith(r"\bibliographystyle"):
            flush_paragraph()
            i += 1
            continue
        if line.startswith(r"\bibliography"):
            flush_paragraph()
            doc.heading("References", 1)
            for key in BIB:
                doc.paragraph(BIB[key].reference_text())
            refs_added = True
            i += 1
            continue
        if line.startswith(r"\section"):
            flush_paragraph()
            m = re.match(r"\\section\*?\{(.+?)\}", line)
            if m:
                doc.heading(latex_to_text(m.group(1)), 1)
            i += 1
            continue
        if line.startswith(r"\subsection"):
            flush_paragraph()
            m = re.match(r"\\subsection\*?\{(.+?)\}", line)
            if m:
                doc.heading(latex_to_text(m.group(1)), 2)
            i += 1
            continue
        if line.startswith(r"\input"):
            flush_paragraph()
            m = re.match(r"\\input\{(.+?)\}", line)
            if m:
                tpath = PROJECT_DIR / m.group(1)
                if tpath.suffix != ".tex":
                    tpath = tpath.with_suffix(".tex")
                caption, rows, note = table_rows_from_tex(tpath)
                doc.paragraph(caption, style="Caption", bold=True)
                doc.add_table(rows)
                if note:
                    doc.paragraph("Note. " + note, style="Caption", italic=True)
            i += 1
            continue
        if line.startswith(r"\begin{figure}"):
            flush_paragraph()
            block = [line]
            i += 1
            while i < len(lines) and not lines[i].strip().startswith(r"\end{figure}"):
                block.append(lines[i])
                i += 1
            if i < len(lines):
                block.append(lines[i])
            process_figure_block(doc, "\n".join(block))
            i += 1
            continue
        if line.startswith("\\") and any(line.startswith(cmd) for cmd in [r"\FloatBarrier", r"\clearpage", r"\appendix", r"\maketitle"]):
            flush_paragraph()
            i += 1
            continue
        paragraph_buf.append(line)
        i += 1
    flush_paragraph()

    if not refs_added:
        doc.heading("References", 1)
        for key in BIB:
            doc.paragraph(BIB[key].reference_text())

    doc.save(OUT_DOCX)


def main() -> None:
    parse_main()
    print(OUT_DOCX)


if __name__ == "__main__":
    main()
