# 🧬 Analisador de Sequências de DNA

[![Abrir no Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/alinesmdados/analisador-dna/blob/main/analisador_dna.ipynb)
![Python](https://img.shields.io/badge/Python-3.9%2B-blue)
![Licença](https://img.shields.io/badge/licen%C3%A7a-MIT-green)

Biblioteca em **Python puro** para análise de sequências de DNA, acompanhada de um notebook que a aplica a dados reais do **NCBI**: o mRNA da insulina humana e o genoma completo do SARS-CoV-2.

Toda a lógica (código genético, busca de ORFs, complemento reverso, etc.) foi implementada do zero, sem Biopython, com o objetivo de demonstrar compreensão dos algoritmos fundamentais da bioinformática.

## Funcionalidades

| Função | O que faz |
|---|---|
| `ler_fasta` / `ler_fasta_texto` | Lê arquivos FASTA com uma ou várias sequências |
| `validar_sequencia` | Verifica se a sequência contém apenas A, C, G, T ou N |
| `contar_nucleotideos` | Contagem de cada base |
| `conteudo_gc` | Percentual de G + C |
| `gc_janela_deslizante` | Perfil de GC ao longo da sequência |
| `assimetria_gc` | GC skew, (G − C)/(G + C), em janela deslizante |
| `complemento_reverso` | Fita complementar no sentido 5'→3' |
| `transcrever` | DNA → RNA |
| `traduzir` | DNA → proteína (código genético padrão) |
| `encontrar_orfs` | ORFs nos 6 quadros de leitura, com coordenadas na fita direta |
| `encontrar_motivo` | Ocorrências de um motivo nas duas fitas (ex.: sítios de restrição) |
| `frequencia_codons` / `frequencia_aminoacidos` | Uso de códons e composição de aminoácidos |
| `temperatura_melting` | Tm estimada de primers (regra de Wallace / fórmula de GC) |
| `resumo` | Relatório com as principais métricas |

## Resultados de validação

O notebook confere a implementação com resultados conhecidos da literatura:

- **Insulina humana (`NM_000207`):** a maior ORF traduzida contém as sequências exatas das cadeias **A** (21 aa) e **B** (30 aa) da insulina madura.
- **SARS-CoV-2 (`NC_045512.2`):** a busca de ORFs recupera a **proteína Spike** completa, com 1.273 aminoácidos, além das demais ORFs do genoma.

## Como usar

### No Google Colab (recomendado)
Clique no botão **"Abrir no Colab"** acima e execute as células em ordem (*Ambiente de execução → Executar tudo*).

### Localmente
```bash
git clone https://github.com/alinesmdados/analisador-dna.git
cd analisador-dna
pip install -r requirements.txt
jupyter notebook analisador_dna.ipynb
```

### Como biblioteca
```python
import analisador_dna as ad

seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
ad.conteudo_gc(seq)              # 56.41
ad.complemento_reverso(seq)      # 'CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT'
ad.traduzir(seq)                 # 'MAIVMGR*KGAR*'
ad.encontrar_orfs(seq, tamanho_minimo_aa=5)
```

## Testes

```bash
pip install pytest
pytest            # testes unitários
python -m doctest analisador_dna.py   # exemplos da documentação
```

## Estrutura

```
analisador-dna/
├── analisador_dna.py         # biblioteca
├── analisador_dna.ipynb      # notebook com as análises
├── test_analisador_dna.py    # testes automatizados
├── requirements.txt
└── LICENSE
```

## Próximos passos

- Alinhamento de sequências (Needleman-Wunsch e Smith-Waterman)
- Suporte a outros códigos genéticos (ex.: mitocondrial)
- Leitura de arquivos FASTQ com escores de qualidade

## Autora

**Aline Machado**: analista de dados, pós-graduada em Ciência de Dados.

## Licença

Distribuído sob a licença MIT. Veja o arquivo [LICENSE](LICENSE).
