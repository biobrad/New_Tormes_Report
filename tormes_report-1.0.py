import os
import pandas as pd
import datapane as dp
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import itertools
from Bio import Phylo
import plotly.figure_factory as ff

# heading stuff

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, '../tormes.log')
def info():
 header=[]
 for word in ['version', 'pipeline']:
  with open(filename, 'r') as fp:
   lines = fp.readlines()
  for line in lines:
   if line.find(word) != -1:
    header.append(line)
 return(header)
turtles = info()
wonk = """# Tormes Report \n""" + "\n" + "### " + turtles[0] + "\n" + "### " + turtles[1] + "\n" + "### " + turtles[2]

## new table creation method

def generate_html(dataframe: pd.DataFrame):
    # get the table HTML from the dataframe
    table_html = dataframe.to_html(table_id="table")
    # construct the complete HTML with jQuery Data tables
    # You can disable paging or enable y scrolling on lines 20 and 21 respectively
    html = f"""
    <html>
    <header>
        <link href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.min.css" rel="stylesheet">
    </header>
    <body>
    {table_html}
    <script src="https://code.jquery.com/jquery-3.6.0.slim.min.js" integrity="sha256-u7e5khyithlIdTpu22PHhENmPcRdFiHRjhAuHcs05RI=" crossorigin="anonymous"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
    <script>
        $(document).ready( function () {{
            $('#table').DataTable({{
                // paging: false,    
                // scrollY: 400,
            }});
        }});
    </script>
    </body>
    </html>
    """
    # return the html
    return html

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
            xanchor="left",
            yanchor="top",
            direction = "left",
            pad={"r": 10, "t": 10},
            x=0.11,
            y=1.16,
            showactive=True,
            buttons=click)])

#generating new html datatable
sarhtml = generate_html(df)

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
taxkrakhtml = generate_html(tax_krak)

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
taxrdphtml = generate_html(tax_rdp)

# MLST 
mlst_data = pd.read_csv("mlst.tab", sep='\t', names=['Sequence', 'Scheme', 'ST', '1', '2', '3', '4', '5', '6', '7'])
mlsthtml = generate_html(mlst_data)

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
    fig.update_layout(width=500, height=700, yaxis={'side': 'right'}, plot_bgcolor='rgba(0,0,0,0)')
    fig.update_traces(hoverinfo='x')
    fig.update_xaxes(ticks="")
    fig.update_yaxes(ticks="", automargin='left+top')
    return(fig)
    
with open ("core_gene_alignment.newick") as core:
    cga = Phylo.read(core, 'newick')
cga = newicktophylo(cga)
with open("accessory_binary_genes.fa.newick") as acc:
    abg = Phylo.read(acc, 'newick')
abg = newicktophylo(abg)

## AMR Resistance Tables
seqs = pd.read_csv("metadata.txt", sep='\t')
farts = seqs['Samples'].tolist()
tablelist=[]
for i in farts:
 res = pd.read_csv(i + '_resfinder.tab', sep='\t')
 res = generate_html(res)
 dfc = pd.read_csv(i + '_card.tab', sep='\t')
 dfc = generate_html(dfc)
 arg = pd.read_csv(i + '_argannot.tab', sep='\t')
 arg = generate_html(arg)
 tablelist.append(dp.Select(label=i, blocks=[dp.Group(dp.HTML(res), label='Resfinder'), dp.Group(dp.HTML(dfc), label='Card'), dp.Group(dp.HTML(arg), label='Argannot')]))
 #tablelist.append(dp.Select(label=i, blocks=[dp.Group(dp.Table(res), label='Resfinder'), dp.Group(dp.Table(dfc), label='Card'), dp.Group(dp.Table(arg), label='Argannot')]))

## AMR Summaries
def makeheatmap(df, color):
   df.set_index('FILE', drop=True, inplace=True)
   df.drop(columns='NUM_FOUND', inplace=True)
   if len(df.columns) < len(df):
    df = df.T
   fig = px.imshow(df, color_continuous_scale=color)
   fig.update_layout(coloraxis_showscale=False, autosize=True)
   fig.update_yaxes(automargin='left+top')
   return(fig)

resfinder = pd.read_csv('resfinder_summary.tab', sep='\t')
resfinder = makeheatmap(resfinder, 'viridis')

card = pd.read_csv('card_summary.tab', sep='\t')
card = makeheatmap(card, 'portland')

argannot = pd.read_csv('argannot_summary.tab', sep='\t')
argannot = makeheatmap(argannot, 'bluered')

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
* Datapane, [Anthias, L *et al., 2022](https://datapane.com/)
"""

report = dp.Report(
    dp.Group(dp.Text("""![](https://github.com/biobrad/New_Tormes_Report/blob/main/tormessmaller.PNG?raw=true)"""), dp.Text(wonk), columns=2),
    dp.Select(type=dp.SelectType.TABS,
        blocks=[
            dp.Group(dp.Text("## Sequencing Assembly Details"), dp.Text(genome_info), dp.Plot(fig), dp.HTML(sarhtml), label="Genome Stats"),
            dp.Group(dp.Text(mlst_info), dp.HTML(mlsthtml), label="MLST"),
            dp.Group(dp.Text(tax_krak_info), dp.HTML(taxkrakhtml), dp.Text(tax_rdp_info), dp.HTML(taxrdphtml),  label="Taxonomy"),
            dp.Group(dp.Text(pan_info), dp.Table(pan_data), dp.Plot(pan_fig), dp.Media(file='./pangenome.png', name='pangenome'), label="Pangenome"),
            dp.Select(label="Phylogenetics", 
                    blocks=[
                        dp.Group(dp.Text('## Core Genes'), dp.Plot(cga), label="Core Genes"),
                        dp.Group(dp.Text('## Accessory Genes'), dp.Plot(abg),  label="Accessory Genes")]),
            dp.Select(label="AMR Summaries",
                    blocks=[dp.Group(dp.Plot(resfinder), label="Resfinder"), 
                            dp.Group(dp.Plot(argannot), label="Argannot"),
                            dp.Group(dp.Plot(card), label="Card")]),
            dp.Select(label="AMR Results",
                    blocks=[*tablelist]),
            dp.Group(dp.Text(citations), label="Citations")
        ]
    ),
)
report.save(path="tormes_report_datapane.html")
