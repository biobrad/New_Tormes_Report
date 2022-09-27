import os
import pandas as pd
import datapane as dp
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import itertools
from Bio import Phylo
import plotly.figure_factory as ff
## Sequence assembly report
df = pd.read_csv("sequencing_assembly_report.txt", sep='\t')
click = []
def butt(val):
    plonk = dict(label=val, method="update", args=[{"y": [df[val]]}])
    click.append(plonk)
head = list(df)
for i in head:
    butt(i)
fig = go.Figure()
fig.add_trace(go.Bar(x=df['ISOLATE'], y=df['GenomeLength'], name="GenomeLength"))
fig.update_layout(title_text='Sequencing Assembly Report', autosize=False, height=500,
    updatemenus=[
        dict(
            type="buttons",
            bgcolor='mediumspringgreen',
            bordercolor='black',
            showactive=True,
            buttons=click)])
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
# MLST 
mlst_data = pd.read_csv("mlst.tab", sep='\t', names=['Sequence', 'Scheme', 'ST', '1', '2', '3', '4', '5', '6', '7'])
mlst_info = """
## Multi-Locus Sequence Typing (MLST)
Further details can be found in [mlst web page](https://github.com/tseemann/mlst).
Symbol | Meaning
-|---------------------------------------------------------------------------------------
~ | novel full length allele similar to match
? | partial match to known allele
`-` | allele missing
"""
# Pangenome
pan_info = """## Pangenome analysis
#### Pangenome genes summary"""
pan_data = pd.read_csv('summary_statistics.txt', sep='\t', names = ['Genes', 'Description', 'Number'])
labels = pan_data['Genes'].tolist()
values = pan_data['Number'].tolist()
pan_fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole =.5)])
#phylogenetics
def newicktophylo(t):
    d = {}
    for x, y in itertools.combinations(t.get_terminals(), 2):
        v = t.distance(x, y)
        d[x.name] = d.get(x.name, {})
        d[x.name][y.name] = v
        d[y.name] = d.get(y.name, {})
        d[y.name][x.name] = v
    for x in t.get_terminals():
        d[x.name][x.name] = 0
    m = pd.DataFrame(d)
    ticknames = m.index.values.tolist()
    fig = ff.create_dendrogram(m, orientation='right', labels=ticknames)
    fig.update_layout(width=800, height=800, yaxis={'side': 'right'}, plot_bgcolor='rgba(0,0,0,0)')
    fig.update_traces(hoverinfo='x')
    fig.update_xaxes(ticks="")
    fig.update_yaxes(ticks="")
    return(fig)
    
with open ("core_gene_alignment.newick") as core:
    cga = Phylo.read(core, 'newick')
cga = newicktophylo(cga)
with open("accessory_binary_genes.fa.newick") as acc:
    abg = Phylo.read(acc, 'newick')
abg = newicktophylo(abg)
## Card Resistance Tables
dirname = '\\Users\\pathe\\Desktop\\Our10Smalto_Plus_NatureStenos\\report_files'
filelist = []
for file in os.listdir(dirname):
    if file.endswith('card.tab'):
        filelist.append(file)
tablelist=[]
for i in filelist:
    seq=i.split('_card')[0]
    df=pd.read_csv(i, sep='\t')
    tablelist.append(dp.Group(dp.DataTable(df), label=seq))

citations = """
### Please cite the following software and databases when using this data for your publication:

* TORMES, [N.M. Quijada *et al*., 2019](https://doi.org/10.1093/bioinformatics/btz220)
* GNU Parallel, [O. Tange, 2018](https://doi.org/10.5281/zenodo.1146014)
* Prinseq, [R. Schmieder and R. Edwards, 2011](https://www.ncbi.nlm.nih.gov/pubmed/21278185)
* SPAdes, [A. Bankevich *et al*., 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342519/)
* QUAST, [A. Gurevich *et al*., 2013](https://www.ncbi.nlm.nih.gov/pubmed/23422339)
* Barrnap, [T. Seemann](https://github.com/tseemann/barrnap)
* Kraken2, [Wood *et al*., 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)
* RDP Classifier, [Q. Wang *et al*., 2007](https://aem.asm.org/content/73/16/5261)
* mlst, [T. Seemann](https://github.com/tseemann/mlst)
* ABRicate, [T. Seemann](https://github.com/tseemann/abricate)
* ResFinder database, [E. Zankari *et al*., 2012](https://www.ncbi.nlm.nih.gov/pubmed/22782487)
* CARD database, [A.G. McArthur *et al*., 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3697360/)
* ARG-ANNOT database, [S.K. Gupta *et al*., 2014](https://www.ncbi.nlm.nih.gov/pubmed/24145532)
* VFDB database, [L. Chen *et al*., 2005](https://www.ncbi.nlm.nih.gov/pubmed/15608208)
* Prodigal, [Hyatt *et al*., 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/)
* Prokka, [T. Seemann, 2014](https://www.ncbi.nlm.nih.gov/pubmed/24642063)
* Roary, [A.J. Page *et al*., 2015](https://www.ncbi.nlm.nih.gov/pubmed/26198102)
* R, [R Development Core Team, 2008](https://www.r-project.org/)
  + ggplot2, [H. Wickham, 2009](https://cran.r-project.org/web/packages/ggplot2/index.html)
  + ggtree, [G. Yu *et al*., 2016](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12628)
  + knitr, [Y. Xie, 2015](https://cran.r-project.org/web/packages/knitr/index.html)
  + plotly, [C. Sievert *et al*., 2017](https://cran.r-project.org/web/packages/plotly/index.html)
  + RColorBrewer, [E. Neuwirth and R.C. Brewer, 2014](https://cran.r-project.org/web/packages/RColorBrewer/index.html)
  + reshape2, [H. Wickham, 2007](https://cran.r-project.org/web/packages/reshape2/index.html)
  + rmarkdown, [J.J. Allaire, 2015](https://cran.r-project.org/web/packages/rmarkdown/index.html)
  + treeio, [L-G. Wang *et al*., 2019](https://academic.oup.com/mbe/article-abstract/37/2/599/5601621?redirectedFrom=fulltext)
"""

report = dp.Report(
    dp.Group(dp.Media(file="./tormessmaller.png", name="Image1"), dp.Text("# Tormes Report \n ### Analysis Performed on 2022-08-26 \n  ### tormes version=1.3.0"), columns=2),
    dp.Select(type=dp.SelectType.TABS,
        blocks=[
            dp.Group(dp.Text("## Sequencing Assembly Details"), dp.Text(genome_info), dp.Plot(fig), dp.DataTable(df), label="Genome Stats"),
            dp.Group(dp.Text(mlst_info), dp.DataTable(mlst_data), label="MLST"),
            dp.Group(dp.Text(tax_krak_info), dp.DataTable(tax_krak), dp.Text(tax_rdp_info), dp.DataTable(tax_rdp),  label="Taxonomy"),
            dp.Group(dp.Text(pan_info), dp.Table(pan_data), dp.Plot(pan_fig), dp.Media(file='./pangenome.png', name='pangenome'), label="Pangenome"),
            dp.Select(label="Phylogenetics", 
                    blocks=[
                        dp.Group(dp.Text('## Core Genes'), dp.Plot(cga), label="Core Genes"),
                        dp.Group(dp.Text('## Accessory Genes'), dp.Plot(abg),  label="Accessory Genes")]
                ),
            dp.Select(label="Card Resistance",
                    blocks=[*tablelist]),
            dp.Group(dp.Text(citations), label="Citations")
        ]
    ),
)
report.save(path="tormes_report1.html")
