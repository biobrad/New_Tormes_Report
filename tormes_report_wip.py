import pandas as pd
import datapane as dp
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

df = pd.read_csv("sequencing_assembly_report.txt", sep='\t')

fig = go.Figure()
fig.add_trace(go.Bar(x=df['ISOLATE'], y=df['GenomeLength'], name="GenomeLength"))

genome_info = """Field | Description
----- | ---------------------------------------------------------------------------------------
**Sample** | Name of the sample
**Reads** | Total number of reads after quality filtering
**AvgReadLen** | Average read length after quality filtering
**Contigs** | Number of contigs of the draft genome (>200bp)
**GenomeLength** | Length (bp) of the draft genome
**LargestContig** | Length (bp) of the largest contig in the genome
**N50** | "Length of the smallest contig in the set that contains the fewest (largest) contigs whose combined length represents at least 50% of the assembly" ([Miller *et al*., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874646/))
**GC** | GC content (%) of the draft genome
**Depth** | Number of times each nucleotide position in the draft genome has a read that align to that position
\* *When already assembled genomes are included for TORMES analysis, "NA" will appear for "Reads", "AvgReadLen" and "Depth" fields.*
"""
tax_krak_info = """### Taxonomic identification by using Kraken2
Taxonomic identification was performed by using Kraken2 ([Wood *et al*., 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)). Further details can be found in the [Kraken2 webpage](https://ccb.jhu.edu/software/kraken2/index.shtml).
Additionally, 16S rRNA genes were extracted from each genome by using [Barrnap](https://github.com/tseemann/barrnap) and used for taxonomic classification by using the RDP Classifier ([Wang *et al*., 2007](https://aem.asm.org/content/73/16/5261)) at a confidence level of 0.8. Further details can be found in the [RDPTools webpage](https://github.com/rdpstaff/classifier).

The number between brackets refers to the percentage of reads (when starting TORMES from raw reads) or contigs (when starting TORMES from already assembled genomes) from each sample covered by the clade rooted at this taxon.
"""

tax_krak = pd.read_csv('taxonomic-identification-kraken2.txt', sep='\t')

tax_rdp_info = """### Taxonomic identification by using RDP Classifier
Field | Description
----- | ---------------------------------------------------------------------------------------
**Sample** | Name of the sample where the 16S rRNA gene was found. Note that the same sample might harbor more than one 16S rRNA gene copy.
**Contig** | Name of the contig where the 16S rRNA gene was found.
**Position** | Position within the contig where the 16S rRNA gene was found (begin-end). Note that the 16S rRNA genes identified for one sample might not harbor the same length (fragmented, *etc.*).
**Strand** | Strand of the 16S rRNA gene defined as forward (+) or reverse (-).
**Order** | Taxonomic order that the 16S rRNA gene was assigned over a confidence level of 0.8.
**Family** | Taxonomic family that the 16S rRNA gene was assigned over a confidence level of 0.8.
**Genus** | Taxonomic genus that the 16S rRNA gene was assigned over a confidence level of 0.8.
"""

tax_rdp = pd.read_csv('taxonomic-identification-16S-rRNA.RDP.txt', sep='\t')

mlst_data = pd.read_csv("mlst.tab", sep='\t')

mlst_info = """
## Multi-Locus Sequence Typing (MLST) {#mlst}

Further details can be found in [mlst web page](https://github.com/tseemann/mlst).

Symbol | Meaning
------ | ---------------------------------------------------------------------------------------
~ | novel full length allele similar to match
? | partial match to known allele
- | allele missing
"""

fig.update_layout(
    updatemenus=[
        dict(
            type="buttons",
            showactive=True,
            buttons=list(
                [
                    dict(
                        label="GenomeLength",
                        method="update",
                        args=[{"y": [df["GenomeLength"]]}],
                    ),
                    dict(
                        label="CONTIGS",
                        method="update",
                        args=[{"y": [df["CONTIGS"]]}],
                    ),
                    dict(
                        label="AvgReadLen",
                        method="update",
                        args=[{"y": [df["AvgReadLen"]]}],
                    ),
                    dict(
                        label="LargestContig",
                        method="update",
                        args=[{"y": [df["LargestContig"]]}],
                    ),
                    dict(
                        label="N50",
                        method="update",
                        args=[{"y": [df["N50"]]}],
                    ),
                    dict(
                        label="GC",
                        method="update",
                        args=[{"y": [df["GC"]]}],
                    ),
                    dict(
                        label="DEPTH",
                        method="update",
                        args=[{"y": [df["DEPTH"]]}],
                    ),
                    dict(
                        label="READS",
                        method="update",
                        args=[{"y": [df["READS"]]}],
                    ),
                ]
            ),
        )
    ]
)
report = dp.Report(
    dp.Group(dp.Media(file="./tormessmaller.png", name="Image1"), dp.Text("# Tormes Report \n ### Analysis Performed on 2022-08-26 \n  ### tormes version=1.3.0"), columns=2),
     dp.Select(
        blocks=[
            dp.Group(dp.Text("## Sequencing Assembly Details"), dp.Text(genome_info), dp.DataTable(df), dp.Plot(fig), label="Genome Stats"),
            dp.Group(dp.Text(mlst_info), dp.DataTable(mlst_data), label="MLST"),
            dp.Group(dp.Text(tax_krak_info), dp.DataTable(tax_krak), dp.Text(tax_rdp_info), dp.DataTable(tax_rdp),  label="Taxonomy"),
            dp.Group()
        ]
    ),
)
report.save(path="tormes_report1.html")
