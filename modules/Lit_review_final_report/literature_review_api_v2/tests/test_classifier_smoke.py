from litreview.models import StudyFinding, Paper
from litreview.classifier import classify_evidence


def finding(**kw):
    base = dict(
        finding_id='f1', assay_type='perturb-seq', target=None, genes=[], pathways=[], phenotypes=[],
        organism='human', tissue_or_cell_type='myeloid cells', disease_or_context=None, summary='test'
    )
    base.update(kw)
    return StudyFinding(**base)


def paper(title, abstract='', pubtype=None):
    return Paper(
        source='pubmed', source_id=title[:10], title=title, abstract=abstract,
        publication_type=pubtype or ['Journal Article']
    )


def test_smoke_cases():
    cases = [
        (finding(pathways=['type I interferon signaling']), paper('Perturb-seq study of type I interferon signaling in human myeloid cells'), 'direct'),
        (finding(pathways=['type I interferon signaling']), paper('Perturb-seq reveals antiviral transcriptional response in human myeloid cells'), 'direct'),
        (finding(pathways=['type I interferon signaling']), paper('Bulk RNA-seq reveals antiviral transcriptional response in human myeloid cells'), 'mechanistic_support'),
        (finding(pathways=['extracellular matrix remodeling']), paper('Perturb-seq analysis of invasion in human myeloid cells'), 'mechanistic_support'),
        (finding(pathways=['type I interferon signaling']), paper('Review of Perturb-seq and type I interferon signaling', pubtype=['Review']), 'background'),
        (finding(pathways=['type I interferon signaling']), paper('Proteomic study of mitochondrial ribosomes in yeast'), 'exclude'),
        (finding(genes=['MET']), paper('Perturb-seq study of metabolic adaptation in human cells'), 'exclude'),
        (finding(genes=['CAT']), paper('Perturb-seq study of catalytic adaptation in human cells'), 'exclude'),
    ]
    for f, p, expected in cases:
        assert classify_evidence(f, p).level.value == expected
