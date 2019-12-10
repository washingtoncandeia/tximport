# Expressão Diferencial em Nível de Transcritos com Kallisto e Sleuth
## Laboratório de Medicina Tropical - IMT - UFRN

Análises utilizando _workflow_ **kallisto -> tximport -> DESeq2**.  

Análises de expressão diferencial de dados de sequenciamento de segunda geração (_high-throughput sequencing reads_, **_RNA-seq_**) realizadas utilizando as ferramentas [kallisto](https://pachterlab.github.io/kallisto/).  

**kallisto** utiliza o conceito de **pseudoalinhamento** para determinar a compatibilidade de _reads_ com seus alvos no genoma, sem que se faça um alinhamento. Quantifica as abundâncias dos transcritos pertencentes a um determinado conjunto de dados de **_RNA-seq_**. Para mais informações, procurar por sua página o artigo do grupo desenvolvedor deste programa:  

Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter, [Near-optimal probabilistic RNA-seq quantification](https://www.nature.com/articles/nbt.3519), Nature Biotechnology **34**, 525–527 (2016), doi:10.1038/nbt.3519. 

Para mais informações acerca dos métodos estatísticos, consultar o artigo do grupo desenvolvedor do programa:  

Harold J. Pimentel, Nicolas Bray, Suzette Puente, Páll Melsted and Lior Pachter, [Differential analysis of RNA-Seq incorporating quantification uncertainty](https://www.nature.com/articles/nmeth.4324), Nature Methods (2017), advanced access http://dx.doi.org/10.1038/nmeth.4324.

-------------------
As análises continuam em andamento no Laboratório de Imunogenética, chefiado pela **Prof. Dr. Selma M.B. Jeronimo**, pertencente ao **Instituto de Medicina Tropical - IMT**, na **Universidade Federal do Rio Grande do Norte**, em Natal.
