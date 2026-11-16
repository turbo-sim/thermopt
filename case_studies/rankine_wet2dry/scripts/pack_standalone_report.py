"""Package the report's eight views, vector figures and audit downloads in one HTML.

Referenced text assets and generated workbooks are embedded once and keyed by
SHA-256. Supplied email and workbook attachments are recorded as citations.
"""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit, urlencode
import base64
import hashlib
import html
import json
import mimetypes
import re

from paths import WEB_ASSETS, OUTPUT


class Attributes(HTMLParser):
    def handle_starttag(self, tag, attrs):
        self.values = dict(attrs)


def attributes(tag):
    parser = Attributes(convert_charrefs=True)
    parser.feed(tag)
    return parser.values


def opening(tag, attrs):
    return "<"+tag+"".join(" "+key if value is None else f' {key}="{html.escape(str(value),quote=True)}"'
                           for key,value in attrs.items())+">"


def download_name(path):
    if path.parent.parent.name == "fluid_comparison":
        return path.parent.name+"_"+path.name
    for index, part in enumerate(path.parts):
        if part.startswith("sensitivity_") and index+1 < len(path.parts):
            prefix = part.removeprefix("sensitivity_")+"_"+path.parts[index+1]
            suffix = {"case.yaml":"input.yaml", "optimal_solution.yaml":"optimum.yaml",
                      "convergence_history.txt":"solver.txt"}.get(path.name)
            if suffix:
                return prefix+"_"+suffix
    return path.name


class Package:
    def __init__(self, pages):
        self.pages = pages
        self.assets = {}
        self.sources = {}
        self.citations = set()

    def asset(self, path):
        path = path.resolve()
        raw = path.read_bytes()
        if path.suffix == ".xlsx":
            encoded = raw
            content = dict(mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                           encoding="base64", data=base64.b64encode(raw).decode("ascii"))
        else:
            text = raw.decode("utf-8")
            if path.suffix == ".svg":
                # Matplotlib's optional external DTD is unnecessary for this image.
                text = re.sub(r"<!DOCTYPE\s+svg\b.*?>", "", text, flags=re.S|re.I)
            mime = {".svg":"image/svg+xml", ".yaml":"application/yaml", ".py":"text/x-python",
                    ".csv":"text/csv", ".json":"application/json", ".txt":"text/plain"}.get(path.suffix, "text/plain")
            encoded = text.encode("utf-8")
            content = dict(mime=mime, text=text)
        key = hashlib.sha256(encoded).hexdigest()
        self.assets.setdefault(key, content)
        self.sources[key] = str(path)
        return key

    def local_path(self, parent, value):
        parsed = urlsplit(html.unescape(value))
        if parsed.scheme or parsed.netloc:
            raise ValueError(f"External resource in standalone report: {value}")
        return (parent / unquote(parsed.path)).resolve(), parsed.fragment

    def body(self, source):
        document = source.read_text(encoding="utf-8")
        match = re.search(r"<body\b[^>]*>(.*)</body>",document,re.S|re.I)
        if not match:
            raise ValueError(f"Missing body in {source}")
        body = match.group(1)

        def image(match):
            attrs = attributes(match.group(0))
            src = attrs.pop("src", None)
            if src:
                path,_ = self.local_path(source.parent,src)
                if path.suffix != ".svg":
                    raise ValueError(f"Expected a vector figure: {path}")
                attrs["data-embedded-image"] = self.asset(path)
            return opening("img",attrs)

        def anchor(match):
            attrs = attributes("<a"+match.group(1)+">")
            content = match.group(2)
            value = attrs.get("href","")
            if not value or value.startswith("#"):
                if value == "#main":
                    attrs["data-focus-main"] = ""
                return opening("a",attrs)+content+"</a>"
            path,fragment = self.local_path(source.parent,value)
            if path in self.pages:
                params = dict(page=path.stem)
                if fragment:
                    params["anchor"] = fragment
                attrs["href"] = "#"+urlencode(params)
                return opening("a",attrs)+content+"</a>"
            if path.suffix == ".eml" or (path.suffix == ".xlsx" and path.name != "optimal_solution.xlsx"):
                self.citations.add(path.name)
                return '<span class="source-citation">'+content+"</span>"
            if path.suffix == ".pdf" and path.name == "report_partial_evaporation.pdf":
                return '<button class="print-action" data-print-page>Print page</button>'
            # SVG remains downloadable. Other figure formats stay in output/;
            # embedding them too would duplicate the same drawing many times.
            if path.suffix in (".png", ".pdf", ".html"):
                return ""
            if path.suffix not in (".svg", ".yaml", ".json", ".csv", ".txt", ".py", ".xlsx"):
                raise ValueError(f"Unsupported embedded download: {path}")
            attrs["href"] = "#download"
            attrs["data-embedded-file"] = self.asset(path)
            attrs["download"] = download_name(path)
            attrs.pop("target",None)
            return opening("a",attrs)+content+"</a>"

        body = re.sub(r"<img\b[^>]*>",image,body,flags=re.I)
        body = re.sub(r"<a\b([^>]*)>(.*?)</a>",anchor,body,flags=re.S|re.I)
        return body


def build(pages_directory, output_file):
    pages_directory,output_file = Path(pages_directory),Path(output_file)
    manifest = json.loads((pages_directory/"pages.json").read_text(encoding="utf-8"))
    page_paths = {(pages_directory/item["path"]).resolve() for item in manifest}
    package = Package(page_paths)
    pages = {}
    for item in manifest:
        path = pages_directory/item["path"]
        pages[path.stem] = dict(title=item["title"],body=package.body(path))
    payload = dict(pages=pages, assets=package.assets)
    data = json.dumps(payload,ensure_ascii=False,separators=(",",":")).replace("</","<\\/")
    css = (WEB_ASSETS/"report.css").read_text(encoding="utf-8")
    css += '\n.print-action{border:0;background:transparent;color:var(--teal);font:inherit;font-weight:600;padding:0;cursor:pointer}.source-citation{font-weight:600}'
    interactions = (WEB_ASSETS/"report.js").read_text(encoding="utf-8")
    router = (WEB_ASSETS/"standalone.js").read_text(encoding="utf-8")
    document = '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>EMPOWER · Partial evaporation</title><style>'+css+'</style><script type="application/json" id="report-payload">'+data+'</script></head><body><div id="report-app"><p style="padding:40px;font:16px Segoe UI">Opening the EMPOWER report…</p></div><noscript>This self-contained report requires JavaScript for navigation and embedded downloads.</noscript><script>window.REPORT_STANDALONE=true;\n'+interactions+'\n'+router+'</script></body></html>'
    output_file.write_text(document,encoding="utf-8")
    packing = dict(pages=len(pages),embedded_assets=len(package.assets),size_bytes=output_file.stat().st_size,
                   original_attachments_embedded=False,source_citations=sorted(package.citations),
                   files={key:source for key,source in package.sources.items()})
    (OUTPUT/"standalone_manifest.json").write_text(json.dumps(packing,indent=2),encoding="utf-8")
    print(f"Standalone report: {len(pages)} views, {len(package.assets)} embedded assets, {output_file.stat().st_size/1e6:.2f} MB")
    return output_file
