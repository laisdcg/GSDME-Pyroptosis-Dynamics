# GSDME-Pyroptosis-Dynamics

## GSDME Systems-Oncology Boolean Network Pipeline

Pipeline em R para integrar o modelo lógico `modelo_GINsim_GSDME_available.zginml` com controle de atratores, perturbações sustentadas, scRNA-seq, perfis multiômicos por paciente e aprendizado por reforço. A versão 1.2.1 inclui leitores e análises reprodutíveis para **GSE125449** e **GSE140228**, concordância observacional ponderada e figuras reorganizadas para publicação. Todas as figuras são produzidas em PDF vetorial e PNG com 600 dpi.

## Escopo

O código preserva as regras e os IDs do modelo GINsim enviado. O modelo contém 27 nós, 39 interações, `GSDME_availability` como único input e `DDR_fixed_ON` como condição mantida em 1 nas análises terapêuticas.

Os cinco módulos metodológicos são:

1. **Controle de atratores e driver nodes:** classifica estruturalmente os nós, enumera intervenções sustentadas de tamanho crescente e mede a proporção de trajetórias que saem do estado resistente e alcançam piroptose sem resistência/proliferação.
2. **Integração single-cell:** normaliza counts, estima para cada gene a probabilidade do estado alto com uma mistura gaussiana de dois componentes e amostra estados iniciais por célula antes da simulação assíncrona.
3. **Validação observacional GEO:** baixa e lê os formatos originais de GSE125449 e GSE140228, mede cobertura dos nós, concordância com endpoints das perturbações, suporte às direções das arestas e destinos Booleanos por grupo celular.
4. **Instâncias multiômicas por paciente:** converte expressão, metilação, CNV e efeitos mutacionais em probabilidades de atividade. Essas probabilidades fazem *soft-clipping* das saídas das regras Booleanas.
5. **Terapia sequencial por aprendizado por reforço:** usa Q-learning tabular para aprender uma sequência temporal orientada especificamente à piroptose mediada por GSDME. A trajetória só termina com `PYROPTOSIS_GSDME = 1`; apoptose sem piroptose não é considerada sucesso terapêutico.

## Avisos científicos essenciais

- Os dados gerados no modo `demo` são **simulados** e servem apenas para verificar o funcionamento do pipeline. Não constituem resultados experimentais.
- “Digital twin” é usado no sentido de **instância computacional personalizada da rede**. O código não produz um gêmeo digital clinicamente validado.
- A busca por driver nodes identifica o menor conjunto **dentro dos nós, direções de intervenção, estados iniciais amostrados e horizonte definidos**. Para uma rede com 27 nós, isso não equivale a uma prova de controle global sobre os 134.217.728 estados possíveis.
- A binarização de scRNA-seq estima atividade transcricional. Expressão de RNA não demonstra automaticamente atividade proteica, clivagem de GSDME, MOMP ou ativação de caspases.
- GSE125449 e GSE140228 não possuem braços experimentais de knockout, overexpression ou tratamento correspondentes às perturbações da rede. Portanto, eles fornecem **concordância observacional**, não confirmação causal das intervenções.
- GSE125449 é analisado nas células anotadas como `Malignant cell`; seus conjuntos 1 e 2 permanecem separados como coortes de descoberta e validação.
- GSE140228 contém células imunes `CD45+`, não células tumorais. Ele é usado somente para suporte do contexto imune/microambiental e nunca como validação direta de uma perturbação intrínseca da célula tumoral.
- miR-204-5p madura, clivagem de GSDME e atividade enzimática de caspases não são medidas de forma confiável por scRNA-seq convencional. O código registra a cobertura real e não imputa esses mecanismos como observados.
- As frequências dos heatmaps são resultados do modelo. A coluna `evidence_scope` diferencia evidência experimental direta em HCC, evidência em outros cânceres e inferência mecanística.

## Requisitos

- R 4.2 ou superior;
- pacotes CRAN: `xml2`, `igraph` e `ggplot2`;
- para os GEOs: `Matrix` e `data.table`;
- conexão à internet na primeira execução, caso os pacotes ainda não estejam instalados.

Instalação manual:

```r
install.packages(c("xml2", "igraph", "ggplot2", "Matrix", "data.table"))
```

## Organização recomendada

Coloque na mesma pasta:

```text
GSDME_systems_oncology_pipeline.R
modelo_GINsim_GSDME_available.zginml
README_GSDME_pipeline.md
```

## Primeira execução

Abra o terminal dentro da pasta e faça primeiro um teste rápido:

```bash
Rscript GSDME_systems_oncology_pipeline.R \
  --model modelo_GINsim_GSDME_available.zginml \
  --out teste_GSDME \
  --quick
```

Execução completa em modo demonstrativo:

```bash
Rscript GSDME_systems_oncology_pipeline.R \
  --model modelo_GINsim_GSDME_available.zginml \
  --out resultados_GSDME \
  --seed 204
```

O script instala os pacotes quando necessário. Para impedir instalação automática, adicione `--no-install`.

## Validação com GSE125449 e GSE140228

Execução recomendada, sem dados demonstrativos:

```bash
Rscript GSDME_systems_oncology_pipeline.R \
  --model modelo_GINsim_GSDME_available.zginml \
  --geo GSE125449,GSE140228 \
  --geo-platform droplet \
  --geo-max-cells 2000 \
  --geo-dir GEO_scRNA_data \
  --out resultados_GEO_GSDME \
  --no-demo \
  --seed 204
```

O download é feito diretamente dos arquivos suplementares oficiais do GEO e reutilizado nas execuções seguintes. O GSE140228 possui uma matriz Droplet comprimida de aproximadamente 250 MB; reserve espaço para a descompressão em memória. `--geo-max-cells` limita as células analisadas depois da leitura da matriz e mantém amostragem estratificada reprodutível.

Para testar o fluxo com menos células:

```bash
Rscript GSDME_systems_oncology_pipeline.R \
  --model modelo_GINsim_GSDME_available.zginml \
  --geo GSE125449,GSE140228 \
  --geo-platform droplet \
  --out teste_GEO_GSDME \
  --no-demo \
  --quick
```

Para analisar a plataforma Smart-seq2 do GSE140228 no lugar da Droplet, use `--geo-platform smartseq2`. Para executar ambas como experimentos separados, use `--geo-platform all`. As plataformas nunca são fundidas antes da normalização, evitando comparação direta indevida de escalas.

### Papel de cada dataset

| Dataset | Conteúdo usado | Papel na análise | Afirmação permitida |
|---|---|---|---|
| GSE125449 | Células `Malignant cell`, Set1 e Set2 separados | Concordância de estados tumorais com endpoints previstos | Suporte observacional em células malignas de câncer hepático |
| GSE140228 | Células imunes CD45+ de tumor, fígado adjacente, sangue, linfonodo e ascite | Contexto imune e heterogeneidade do microambiente | Suporte microambiental; não validação tumoral direta |

Para cada experimento, o código calcula:

1. taxa de detecção de cada nó mapeável da rede;
2. probabilidade média de ativação transcricional;
3. distância absoluta ponderada entre o perfil observado e cada endpoint simulado de perturbação, convertida em escore de concordância entre 0 e 1;
4. correlação de Spearman para arestas cujos nós de origem e destino são observáveis;
5. probabilidades de destinos Booleanos por grupo celular.

O escore principal de concordância é:

```text
concordância_ponderada =
  1 − soma[peso × |probabilidade_observada − frequência_simulada|] / soma(pesos)
```

Os pesos são 4 para o nó diretamente perturbado, 2 para seus alvos regulatórios de primeira ordem e 1 para os demais nós mensuráveis. O script também preserva o escore não ponderado, registra quais alvos diretos foram observados e classifica como `Indirect only` os casos em que o alvo perturbado não foi medido. O número e a identidade dos nós comparados são informados nas tabelas. Um escore alto indica semelhança de estado, não efeito causal do tratamento.

## Uso com scRNA-seq real

O arquivo pode ser CSV, TSV/TXT ou RDS. A matriz deve conter genes nas linhas e células nas colunas.

Exemplo CSV:

```text
gene,Cell_001,Cell_002,Cell_003
GSDME,10,0,4
CASP3,8,5,7
SIRT1,0,12,3
MALAT1,1,15,4
```

Execute:

```bash
Rscript GSDME_systems_oncology_pipeline.R \
  --model modelo_GINsim_GSDME_available.zginml \
  --scrna counts_scRNA.csv \
  --out resultados_scRNA \
  --no-demo
```

Aliases reconhecidos incluem:

| Nó lógico | Genes/aliases aceitos |
|---|---|
| `GSDME_availability` | GSDME, DFNA5 |
| `CYTOCHROME_C` | CYCS, CYTOCHROME_C |
| `PUMA` | BBC3, PUMA |
| `p53_ACTIVE` | TP53, P53 |
| `miR_204_5p` | MIR204, MIR204-5P, MIR_204_5P |
| `PGC1A` | PPARGC1A, PGC1A |
| `p21` | CDKN1A, P21 |
| `CyclinD_CDK46` | CCND1, CDK4, CDK6 |

Quando mais de um gene é associado ao mesmo nó, o código usa a maior probabilidade de atividade. Essa regra é configurável no objeto `gene_map`.

## Uso com multiômicas reais

Use uma tabela longa, com uma linha por combinação paciente–nó:

```text
patient_id,node,expression_z,methylation_beta,cnv_log2,mutation_effect
P01,GSDME_availability,-1.2,0.85,-0.30,0
P01,SIRT1,1.5,0.20,0.40,0
P01,p53_ACTIVE,-0.8,0.60,-0.20,-1
P02,GSDME_availability,1.1,0.15,0.20,0
```

Definições:

- `expression_z`: expressão padronizada dentro da coorte;
- `methylation_beta`: valor entre 0 e 1;
- `cnv_log2`: razão de cópia em escala log2;
- `mutation_effect`: `-1` para perda de função, `0` para ausência/efeito desconhecido e `1` para ganho de função.

Execute:

```bash
Rscript GSDME_systems_oncology_pipeline.R \
  --model modelo_GINsim_GSDME_available.zginml \
  --scrna counts_scRNA.csv \
  --multiomics pacientes_multiomics.csv \
  --out resultados_reais \
  --no-demo
```

O prior multiômico é calculado por:

```text
logit(P[nó ativo]) =
  0.9 × expressão_z
  − 2.0 × (metilação_beta − 0.5)
  + 1.1 × cnv_log2
  + 4.0 × mutation_effect
```

Durante a atualização:

```text
P(final = 1) =
  (1 − λ) × saída_da_regra_Booleana
  + λ × prior_multiômico
```

O valor padrão é `λ = 0.35`. Esses pesos são hipóteses configuráveis e devem ser calibrados com uma coorte de treinamento independente antes de qualquer interpretação preditiva.

## Perturbações incluídas

Todas as perturbações são simuladas sob `DDR_fixed_ON = 1`.

| Perturbação | Interpretação |
|---|---|
| GSDME KO | Testa a mudança de piroptose para apoptose quando CASP3 está ativa |
| CASP3 KO | Testa a dependência da clivagem de GSDME por caspase-3 |
| miR-204-5p OE/KO | Testa a modulação de SIRT1 |
| MALAT1 OE/KO | Testa o sequestro de miR-204 e liberação de SIRT1 |
| SIRT1 OE/KO | Testa sobrevivência/resistência versus morte mitocondrial |
| BCL2 KO | Testa a liberação de BAX |

O arquivo `04_scientific_perturbation_evidence.csv` registra para cada intervenção o resultado esperado, o escopo da evidência e PMID/DOI. A direção inversa de uma interação validada não é automaticamente uma validação experimental direta; por isso essas condições aparecem como “mechanistic inverse inferred”.

## Saídas

### Figuras principais

Cada figura é gerada em PDF vetorial e PNG de 600 dpi:

1. `Figure_01_GSDME_logical_network`: rede com ativação em verde, inibição em vermelho, input e DDR fixo destacados;
2. `Figure_02_in_silico_perturbation_heatmap`: frequência de ativação dos fenótipos nas perturbações curadas, explicitamente identificada como resultado do modelo;
3. `Figure_03_minimum_driver_nodes`: melhores intervenções candidatas, identificadas por nome biológico e estado ON/OFF, sem legenda redundante de tamanho igual a 1;
4. `Figure_04_GEO_node_detection_heatmap`: cobertura dos nós em cada experimento GEO;
5. `Figure_05_GEO_weighted_perturbation_concordance`: concordância observacional ponderada com os endpoints simulados;
6. `Figure_06_GEO_projected_fate_composition`: potenciais de destino projetados pelo modelo por coorte ou grupo imune;
7. `Figure_07_pyroptosis_oriented_RL_therapy`: trajetória da política aprendida por Q-learning com alvo específico em piroptose mediada por GSDME.

Assim, uma execução GEO com `--no-demo` produz a sequência principal completa de **Figure 01 a Figure 07**, sem lacunas entre os números.

### Figuras suplementares condicionais

- `Supplementary_Figure_S1_GEO_regulatory_edge_support`: criada somente quando existem pelo menos três arestas com correlação avaliável; quando a cobertura é insuficiente, a tabela é mantida e a figura esparsa é omitida;
- `Supplementary_Figure_S2_generic_scRNA_cell_fates`: probabilidades por célula para uma matriz fornecida com `--scrna` ou no modo demonstrativo;
- `Supplementary_Figure_S3_generic_scRNA_composition`: composição de destinos para a matriz single-cell genérica;
- `Supplementary_Figure_S4_multiomic_digital_twins`: probabilidades fenotípicas por paciente quando `--multiomics` é fornecido ou no modo demonstrativo.

As figuras usam paleta distinguível para daltonismo, rótulos biológicos legíveis e dimensões ajustadas para reduzir espaço vazio. Os PDFs são preferíveis para montagem final do artigo porque permanecem vetoriais.

### Tabelas

- nós, regras e arestas importados do GINsim;
- atratores de referência resistente e piroptótico;
- evidências das perturbações;
- frequências de ativação por perturbação;
- ranking estrutural e conjuntos mínimos candidatos;
- probabilidades single-cell;
- priors e fenótipos multiômicos por paciente;
- sequência aprendida pelo agente;
- manifesto dos arquivos GEO baixados, escopo interpretativo e afirmações permitidas;
- cobertura de nós, concordância ponderada e não ponderada, disponibilidade do alvo direto, grau de evidência e probabilidades de destino nos GEOs;
- suporte correlacional às arestas regulatórias e resumo separado por dataset;
- resumo do aprendizado por reforço informando explicitamente se o alvo piroptótico foi alcançado dentro do horizonte;
- manifesto dos arquivos.

O diretório `logs` contém o registro de execução e `sessionInfo.txt`.

## Solução de problemas

Se uma versão anterior mostrar vários inputs como `NA` e interromper em `list2env`, substitua o script pela versão atual 1.2.1. Nos arquivos GINsim, o atributo `input` normalmente está ausente nos nós internos; o importador corrigido interpreta apenas `input="true"` como input e preserva explicitamente os nomes dos nós ao criar o ambiente das regras lógicas.

Para não misturar imagens antigas com a nova numeração, use um novo diretório em `--out`, por exemplo `resultados_GEO_GSDME_v120`.

## Interpretação do módulo de controle

O estado maligno de referência é `SURVIVAL = 1` e `RESISTANCE = 1` sob DDR ativo. O alvo padrão é:

```text
PYROPTOSIS_GSDME = 1
RESISTANCE = 0
PROLIFERATION = 0
```

A busca testa intervenções sustentadas em MALAT1, miR-204-5p, SIRT1, PGC-1α, BCL2, BAX, p53 e CASP3. Ela começa com intervenções únicas e aumenta o tamanho até encontrar conjuntos com sucesso mínimo de 95% nas trajetórias amostradas.

Para um artigo, descreva o resultado como **robust candidate control set under the tested state and intervention space**, e não como garantia universal, a menos que posteriormente seja aplicada verificação simbólica/exaustiva.

## Interpretação do aprendizado por reforço

As ações disponíveis são:

```text
NONE
MALAT1_KO
miR204_OE
SIRT1_KO
PGC1A_KO
BCL2_KO
```

O agente recebe a maior recompensa somente por `PYROPTOSIS_GSDME = 1`, recompensas intermediárias menores por GSDME-N, CASP3 e BAX, e penalidades por resistência, proliferação, sobrevivência, apoptose sem piroptose e custo de intervenção. O arquivo `20_RL_pyroptosis_target_summary.csv` declara se o alvo foi alcançado. Se ele não for atingido dentro do horizonte, a própria legenda da Figure 07 informa que não há demonstração de sucesso. A sequência é uma hipótese computacional, não uma recomendação terapêutica.

## Referências metodológicas e biológicas

1. Jiang G, Wen L, Zheng H, Jian Z, Deng W. *miR-204-5p targeting SIRT1 regulates hepatocellular carcinoma progression*. Cell Biochem Funct. 2016;34:505–510. [PMID 27748572](https://pubmed.ncbi.nlm.nih.gov/27748572/) — DOI: 10.1002/cbf.3223.
2. Hou Z et al. *The long non-coding RNA MALAT1 promotes the migration and invasion of hepatocellular carcinoma by sponging miR-204 and releasing SIRT1*. Tumour Biol. 2017. [PMID 28720061](https://pubmed.ncbi.nlm.nih.gov/28720061/) — DOI: 10.1177/1010428317718135.
3. Wang Y et al. *Chemotherapy drugs induce pyroptosis through caspase-3 cleavage of a gasdermin*. Nature. 2017;547:99–103. [PMID 28459430](https://pubmed.ncbi.nlm.nih.gov/28459430/) — DOI: 10.1038/nature22393.
4. Sun X et al. *Germacrone induces caspase-3/GSDME activation and enhances ROS production, causing HepG2 pyroptosis*. Exp Ther Med. 2022;24:456. [PMID 35747157](https://pubmed.ncbi.nlm.nih.gov/35747157/) — DOI: 10.3892/etm.2022.11383.
5. Hou W, Tamura T, Ching WK, Akutsu T. *Finding and analyzing the minimum set of driver nodes in control of Boolean networks*. Advances in Complex Systems. 2016;19:1650006. DOI: 10.1142/S0219525916500065.
6. Yang G et al. *Target Control in Logical Models Using the Domain of Influence of Nodes*. Front Physiol. 2018;9:454. [Article](https://www.frontiersin.org/journals/physiology/articles/10.3389/fphys.2018.00454/full).
7. Magaña-López G et al. *scBoolSeq: Linking scRNA-seq statistics and Boolean dynamics*. PLoS Comput Biol. 2024;20:e1011620. [Article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011620).
8. Papagiannis G et al. *Deep Reinforcement Learning for Control of Probabilistic Boolean Networks*. 2021. [Repository record](https://openresearch.surrey.ac.uk/esploro/outputs/conferencePaper/Deep-Reinforcement-Learning-for-Control-of/99522222202346).
9. Ma L, Hernandez MO, et al. *Tumor Cell Biodiversity Drives Microenvironmental Reprogramming in Liver Cancer*. Cancer Cell. 2019;36:418–430.e6. [PMID 31588021](https://pubmed.ncbi.nlm.nih.gov/31588021/) — DOI: 10.1016/j.ccell.2019.08.007; GEO: GSE125449.
10. Zhang Q et al. *Landscape and Dynamics of Single Immune Cells in Hepatocellular Carcinoma*. Cell. 2019;179:829–845.e20. [PMID 31675496](https://pubmed.ncbi.nlm.nih.gov/31675496/) — DOI: 10.1016/j.cell.2019.10.003; GEO: GSE140228.

## Reprodutibilidade

- A semente padrão é `204` e pode ser alterada com `--seed`.
- Não edite os IDs no GINsim sem atualizar as regras e o `gene_map`.
- Guarde o `.zginml`, o script, o README, as tabelas de entrada e `sessionInfo.txt` junto ao material suplementar.
- Preserve os arquivos GEO em `GEO_scRNA_data` para que a análise possa ser repetida sem novo download.
- Não agregue diretamente GSE125449 e GSE140228: eles representam compartimentos celulares e tecnologias diferentes. Compare apenas métricas resumidas e seus respectivos escopos.
- Para resultados finais, aumente repetições somente depois de validar o modo `--quick`.
