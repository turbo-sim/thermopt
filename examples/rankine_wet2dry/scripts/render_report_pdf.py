"""Render report_content.json using reportlab (no thermodynamic dependencies)."""
from pathlib import Path
import json
from xml.sax.saxutils import escape
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, KeepTogether
from PIL import Image as PILImage

from paths import STUDY as HERE, RESULTS


def main():
    fontdir = Path("C:/Windows/Fonts")
    pdfmetrics.registerFont(TTFont("Report", str(fontdir / "segoeui.ttf")))
    pdfmetrics.registerFont(TTFont("ReportBold", str(fontdir / "segoeuib.ttf")))
    pdfmetrics.registerFontFamily("Report", normal="Report", bold="ReportBold", italic="Report", boldItalic="ReportBold")
    styles = getSampleStyleSheet()
    for key in ("Normal", "BodyText", "Heading1", "Heading2", "Heading3", "Title"):
        styles[key].fontName = "ReportBold" if key.startswith("Heading") or key == "Title" else "Report"
    styles["BodyText"].fontSize = 9
    styles["BodyText"].leading = 13.5
    styles["BodyText"].spaceAfter = 8
    styles["Heading2"].spaceBefore = 18
    styles.add(ParagraphStyle("Cell", fontName="Report", fontSize=7, leading=9))
    styles.add(ParagraphStyle("HeaderCell", fontName="ReportBold", fontSize=7, leading=9, textColor=colors.white))
    styles.add(ParagraphStyle("CaptionText", fontName="Report", fontSize=8, leading=11, textColor=colors.HexColor("#526777"), spaceAfter=10))
    nodes = json.loads((RESULTS / "report_content.json").read_text(encoding="utf-8"))
    story = []
    width = 174*mm
    skip = False
    for index, node in enumerate(nodes):
        if skip:
            skip = False
            continue
        kind = node["kind"]
        if kind == "heading":
            title = Paragraph(escape(node["text"]), styles["Title"] if node["level"] == 1 else styles["Heading2"])
            if index+1 < len(nodes) and nodes[index+1]["kind"] == "paragraph":
                lead = Paragraph(escape(nodes[index+1]["text"]), styles["BodyText"])
                story.append(KeepTogether([title, lead]))
                skip = True
            else:
                story.append(title)
        elif kind == "paragraph":
            story.append(Paragraph(escape(node["text"]), styles["BodyText"]))
        elif kind == "table":
            rows = [[Paragraph(escape(str(v)), styles["HeaderCell"] if i == 0 else styles["Cell"]) for v in row] for i,row in enumerate([node["headers"]]+node["rows"])]
            n = len(node["headers"])
            weights = [2]+[1]*(n-1) if n >= 5 else [1.2]+[2]*(n-1) if n == 2 else [1.6]+[1]*(n-1)
            t = Table(rows, colWidths=[width*w/sum(weights) for w in weights], repeatRows=1, hAlign="LEFT")
            t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),colors.HexColor("#163e58")), ("VALIGN",(0,0),(-1,-1),"TOP"), ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white, colors.HexColor("#f3f7fa")]), ("LINEBELOW",(0,0),(-1,-1),.3,colors.HexColor("#d8e2e7")), ("TOPPADDING",(0,0),(-1,-1),5), ("BOTTOMPADDING",(0,0),(-1,-1),5)]))
            story.extend([t, Spacer(1, 8)])
        elif kind == "figure":
            path = HERE / (node["path"]+".png")
            with PILImage.open(path) as im:
                w,h = im.size
            scale = min(width/w, 212*mm/h)
            story.append(KeepTogether([Image(str(path), width=w*scale, height=h*scale), Paragraph(escape(node["caption"]),styles["CaptionText"])]))
    def footer(canvas, doc):
        canvas.setFont("Report", 8)
        canvas.setFillColor(colors.HexColor("#526777"))
        canvas.drawString(18*mm, 12*mm, "EMPOWER | Partial evaporation | Updated boundary conditions")
        canvas.drawRightString(192*mm, 12*mm, str(doc.page))
    output = HERE / "report_partial_evaporation.pdf"
    doc = SimpleDocTemplate(str(output), pagesize=(210*mm,297*mm), rightMargin=18*mm,leftMargin=18*mm,topMargin=17*mm,bottomMargin=20*mm,
                            title="Partial evaporation in geothermal Rankine cycles",author="ThermoOpt")
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(output)


if __name__ == "__main__":
    main()
