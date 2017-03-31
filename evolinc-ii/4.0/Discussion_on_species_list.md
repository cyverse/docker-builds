# Best practices for preparing species list file:

To perform its phylogenetic comparisons, Evolinc-II needs a list of all the species you are searching through, in order of phylogenetic relationship to the query species. This list uses the same four letter Genus/Species abbreviation that the BLASTing list does (first letter of genus followed by first three letters of species). For instance, humans would be **Hsap**, mouse would be **Mmus**, *Arabidopsis thaliana* would be **Atha**, etc. When setting up the species list file, it should be a single column file that starts with the query species, using the four letter code. Let's start by using the the mustard (Brassicaceae) lineage as an example, particularly the ten species listed in the Evolinc manuscript [preprint](http://biorxiv.org/content/early/2017/02/20/110148) (*Arabidopsis thaliana*, *Arabidopsis lyrata*, *Capsella rubella*, *Leavenworthia alabamica*, *Brassica rapa*, *Brassica oleraceae*, *Schrenkiella parvula*, *Eutrema salsugineu*, *Aethionema arabicum*, and *Tarenaya hassleriana*). I have listed them out in order of relationship to *Arabidopsis thaliana*. If you don't know the relationship of the species you are examining, you can usually google "Species X phylogeny" or search for them on [wikipedia](https://www.wikipedia.org) to get a starting point. Here is the resulting species list for these ten species:

```
Atha
Alyr
Crub
Lala
Brap
Bole
Spar
Esal
Aara
Thas
```

If you have species that diverged from your query at the same node, then you can list them in any order (because they are all equidistant to the query). In the above example, **Brap**, **Bole**, **Spar**, and **Esal** all diverged from **Atha** at the same node. Thus, the following list would be equivalent:

```
Atha
Alyr
Crub
Lala
Esal
Spar
Brap
Bole
Aara
Thas
```

If you are running the RAxML portion of Evolinc-II, then you will have generated a newick formatted species-tree (using the four-letter abbreviations) describing the relationship of all species to one another. The order that you list species in the species-tree file (starting with the query) is an acceptable way to write out the species list column.

For those more comfortable with a mammalian lineage, here is the mammalian species list to go along with the TERC analysis shown in the Evolinc manuscript [preprint](http://biorxiv.org/content/early/2017/02/20/110148):

```
Hsap
Ptro
Ggor
Pabe
Nleu
Mmul
Ogar
Mmur
Efla
Ocun
Cpor
Itri
Mmus
Rnor
Mluc
Ecab
Fcat
Mput
Clup
Oari
Sscr
Hgla
Lafr
Tman
Mdom
Oana
```