"""Locate arXiv preprints for three paywalled references by title search."""
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import io
import os

UA = {'User-Agent': 'PaperA-ref-audit/1.0 (mailto:8073592@qq.com)'}
NS = {'a': 'http://www.w3.org/2005/Atom',
      'arxiv': 'http://arxiv.org/schemas/atom'}

QUERIES = [
    'ti:"Ollivier curvature of random geometric graphs"',
    'ti:"Ricci curvature and eigenvalue estimate on locally finite graphs"',
    'ti:"Ollivier-Ricci idleness functions"',
]

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "ref_verify", "pdfs")


def main():
    for q in QUERIES:
        url = ("http://export.arxiv.org/api/query?search_query="
               + urllib.parse.quote(q) + "&max_results=3")
        data = urllib.request.urlopen(
            urllib.request.Request(url, headers=UA), timeout=60).read()
        root = ET.fromstring(data.decode('utf-8', 'replace'))
        print('###', q)
        for e in root.findall('a:entry', NS):
            title = ' '.join((e.find('a:title', NS).text or '').split())
            aid = e.find('a:id', NS).text.split('/abs/')[-1]
            jr = e.find('arxiv:journal_ref', NS)
            ad = e.find('arxiv:doi', NS)
            au = [(x.text or '') for x in e.findall('a:author/a:name', NS)]
            print('   id:', aid, '|', title[:72])
            print('      ', ', '.join(au))
            print('       jref:', (jr.text if jr is not None else None),
                  '| doi:', (ad.text if ad is not None else None))
            # download the first hit per query
            pdf = os.path.join(OUT, aid.replace('/', '_') + '.pdf')
            try:
                d = urllib.request.urlopen(urllib.request.Request(
                    'https://arxiv.org/pdf/' + aid, headers=UA),
                    timeout=90).read()
                io.open(pdf, 'wb').write(d)
                print('       downloaded ->', os.path.basename(pdf), len(d))
            except Exception as exc:
                print('       dl err', exc)
            break
        print()


if __name__ == "__main__":
    main()
