"""
analisador_dna.py
=================

Biblioteca em Python puro para análise básica de sequências de DNA.

Funcionalidades:
    - Leitura de arquivos FASTA (uma ou várias sequências)
    - Validação e limpeza de sequências
    - Contagem de nucleotídeos e conteúdo GC (global e em janela deslizante)
    - Complemento reverso
    - Transcrição (DNA -> RNA) e tradução (DNA -> proteína)
    - Busca de ORFs (Open Reading Frames) nos 6 quadros de leitura
    - Busca de motivos (subsequências) em ambas as fitas
    - Frequência de códons
    - Temperatura de melting (Tm) para primers

Não depende de bibliotecas externas: toda a lógica foi implementada
do zero, com fins didáticos.

Autora: Aline Machado
Licença: MIT
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------

NUCLEOTIDEOS_VALIDOS = set("ACGTN")

COMPLEMENTO = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Código genético padrão (NCBI, tabela 1). '*' representa códon de parada.
_BASES = "TCAG"
_AMINOACIDOS = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
CODIGO_GENETICO = {
    a + b + c: _AMINOACIDOS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}

CODON_INICIO = "ATG"
CODONS_PARADA = {"TAA", "TAG", "TGA"}


# ---------------------------------------------------------------------------
# Leitura e validação
# ---------------------------------------------------------------------------

@dataclass
class RegistroFasta:
    """Uma sequência lida de um arquivo FASTA."""

    identificador: str
    descricao: str
    sequencia: str

    def __len__(self) -> int:
        return len(self.sequencia)


def ler_fasta_texto(texto: str) -> list[RegistroFasta]:
    """Lê sequências a partir de um texto no formato FASTA.

    >>> regs = ler_fasta_texto(">seq1 exemplo\\nATGC\\nGGCC\\n")
    >>> regs[0].identificador, regs[0].sequencia
    ('seq1', 'ATGCGGCC')
    """
    registros: list[RegistroFasta] = []
    cabecalho = None
    partes: list[str] = []

    for linha in texto.splitlines():
        linha = linha.strip()
        if not linha:
            continue
        if linha.startswith(">"):
            if cabecalho is not None:
                registros.append(_montar_registro(cabecalho, partes))
            cabecalho = linha[1:].strip()
            partes = []
        else:
            if cabecalho is None:
                raise ValueError("Arquivo FASTA inválido: sequência antes do cabeçalho '>'.")
            partes.append(linha)

    if cabecalho is not None:
        registros.append(_montar_registro(cabecalho, partes))

    return registros


def _montar_registro(cabecalho: str, partes: list[str]) -> RegistroFasta:
    identificador, _, descricao = cabecalho.partition(" ")
    return RegistroFasta(identificador, descricao, limpar_sequencia("".join(partes)))


def ler_fasta(caminho: str) -> list[RegistroFasta]:
    """Lê um arquivo FASTA do disco."""
    with open(caminho, encoding="utf-8") as arquivo:
        return ler_fasta_texto(arquivo.read())


def limpar_sequencia(seq: str) -> str:
    """Remove espaços, números e quebras de linha e converte para maiúsculas."""
    return "".join(c for c in seq.upper() if c.isalpha())


def validar_sequencia(seq: str) -> bool:
    """Retorna True se a sequência contém apenas A, C, G, T ou N."""
    return bool(seq) and set(seq.upper()) <= NUCLEOTIDEOS_VALIDOS


def _exigir_valida(seq: str) -> str:
    seq = limpar_sequencia(seq)
    if not validar_sequencia(seq):
        invalidos = sorted(set(seq) - NUCLEOTIDEOS_VALIDOS)
        raise ValueError(f"Sequência de DNA inválida. Caracteres não permitidos: {invalidos}")
    return seq


# ---------------------------------------------------------------------------
# Composição
# ---------------------------------------------------------------------------

def contar_nucleotideos(seq: str) -> dict[str, int]:
    """Conta A, C, G, T (e N, se houver).

    >>> contar_nucleotideos("AACGTT")
    {'A': 2, 'C': 1, 'G': 1, 'T': 2}
    """
    seq = _exigir_valida(seq)
    contagem = Counter(seq)
    resultado = {base: contagem.get(base, 0) for base in "ACGT"}
    if contagem.get("N"):
        resultado["N"] = contagem["N"]
    return resultado


def conteudo_gc(seq: str) -> float:
    """Percentual de G + C na sequência (ignora N no denominador).

    >>> conteudo_gc("GGCCAATT")
    50.0
    """
    seq = _exigir_valida(seq)
    validas = len(seq) - seq.count("N")
    if validas == 0:
        return 0.0
    return round(100 * (seq.count("G") + seq.count("C")) / validas, 2)


def gc_janela_deslizante(seq: str, janela: int = 100, passo: int = 10) -> list[tuple[int, float]]:
    """Calcula o conteúdo GC ao longo da sequência.

    Retorna uma lista de tuplas (posição central da janela, %GC).
    Útil para identificar ilhas CpG e regiões ricas em GC.
    """
    seq = _exigir_valida(seq)
    if janela <= 0 or passo <= 0:
        raise ValueError("janela e passo devem ser positivos.")
    if len(seq) < janela:
        return [(len(seq) // 2, conteudo_gc(seq))]
    return [
        (inicio + janela // 2, conteudo_gc(seq[inicio:inicio + janela]))
        for inicio in range(0, len(seq) - janela + 1, passo)
    ]


def assimetria_gc(seq: str, janela: int = 1000, passo: int = 100) -> list[tuple[int, float]]:
    """GC skew = (G - C) / (G + C) em janela deslizante.

    Em genomas bacterianos, a mudança de sinal do GC skew costuma
    indicar a origem e o término de replicação.
    """
    seq = _exigir_valida(seq)
    resultado = []
    for inicio in range(0, max(len(seq) - janela + 1, 1), passo):
        trecho = seq[inicio:inicio + janela]
        g, c = trecho.count("G"), trecho.count("C")
        skew = (g - c) / (g + c) if (g + c) else 0.0
        resultado.append((inicio + len(trecho) // 2, round(skew, 4)))
    return resultado


# ---------------------------------------------------------------------------
# Transformações
# ---------------------------------------------------------------------------

def complemento_reverso(seq: str) -> str:
    """Retorna a fita complementar lida de 5' para 3'.

    >>> complemento_reverso("ATGCCG")
    'CGGCAT'
    """
    return _exigir_valida(seq).translate(COMPLEMENTO)[::-1]


def transcrever(seq: str) -> str:
    """Transcreve DNA (fita codificadora) em RNA mensageiro: T -> U.

    >>> transcrever("ATGCTT")
    'AUGCUU'
    """
    return _exigir_valida(seq).replace("T", "U")


def traduzir(seq: str, parar_no_stop: bool = False) -> str:
    """Traduz DNA em sequência de aminoácidos (código de uma letra).

    Códons com N viram 'X'. Códons de parada viram '*'.

    >>> traduzir("ATGGCCTAA")
    'MA*'
    >>> traduzir("ATGGCCTAAGGG", parar_no_stop=True)
    'MA'
    """
    seq = _exigir_valida(seq)
    proteina = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODIGO_GENETICO.get(seq[i:i + 3], "X")
        if aa == "*" and parar_no_stop:
            break
        proteina.append(aa)
    return "".join(proteina)


# ---------------------------------------------------------------------------
# ORFs
# ---------------------------------------------------------------------------

@dataclass
class ORF:
    """Uma Open Reading Frame encontrada na sequência."""

    fita: str          # '+' (direta) ou '-' (reversa)
    quadro: int        # 1, 2 ou 3
    inicio: int        # posição inicial (1-based, na fita direta)
    fim: int           # posição final (1-based, inclusiva, na fita direta)
    proteina: str

    @property
    def tamanho_nt(self) -> int:
        return self.fim - self.inicio + 1

    @property
    def tamanho_aa(self) -> int:
        return len(self.proteina)


def encontrar_orfs(seq: str, tamanho_minimo_aa: int = 100) -> list[ORF]:
    """Busca ORFs (ATG ... códon de parada) nos 6 quadros de leitura.

    Para cada códon de parada, considera apenas o ATG mais a montante
    (a ORF mais longa), como faz o ORFfinder do NCBI.

    As coordenadas retornadas são sempre relativas à fita direta (1-based).
    """
    seq = _exigir_valida(seq)
    n = len(seq)
    orfs: list[ORF] = []

    for fita, s in (("+", seq), ("-", complemento_reverso(seq))):
        for quadro in range(3):
            inicio_atual = None
            for i in range(quadro, n - 2, 3):
                codon = s[i:i + 3]
                if inicio_atual is None and codon == CODON_INICIO:
                    inicio_atual = i
                elif inicio_atual is not None and codon in CODONS_PARADA:
                    proteina = traduzir(s[inicio_atual:i])
                    if len(proteina) >= tamanho_minimo_aa:
                        fim_local = i + 2  # inclui o códon de parada
                        if fita == "+":
                            ini, fim = inicio_atual + 1, fim_local + 1
                        else:
                            ini, fim = n - fim_local, n - inicio_atual
                        orfs.append(ORF(fita, quadro + 1, ini, fim, proteina))
                    inicio_atual = None

    return sorted(orfs, key=lambda o: o.tamanho_aa, reverse=True)


# ---------------------------------------------------------------------------
# Motivos e códons
# ---------------------------------------------------------------------------

def encontrar_motivo(seq: str, motivo: str, ambas_fitas: bool = True) -> list[tuple[int, str]]:
    """Encontra todas as ocorrências (inclusive sobrepostas) de um motivo.

    Retorna lista de (posição 1-based na fita direta, fita).

    >>> encontrar_motivo("GAATTCGAATTC", "GAATTC", ambas_fitas=False)
    [(1, '+'), (7, '+')]
    """
    seq = _exigir_valida(seq)
    motivo = _exigir_valida(motivo)
    achados = _ocorrencias(seq, motivo, "+")

    if ambas_fitas:
        motivo_rc = complemento_reverso(motivo)
        if motivo_rc != motivo:  # evita duplicar motivos palindrômicos
            achados += _ocorrencias(seq, motivo_rc, "-")

    return sorted(achados)


def _ocorrencias(seq: str, motivo: str, fita: str) -> list[tuple[int, str]]:
    posicoes = []
    inicio = seq.find(motivo)
    while inicio != -1:
        posicoes.append((inicio + 1, fita))
        inicio = seq.find(motivo, inicio + 1)
    return posicoes


def frequencia_codons(seq: str) -> dict[str, int]:
    """Conta os códons de uma sequência codificadora (quadro 1)."""
    seq = _exigir_valida(seq)
    return dict(Counter(
        seq[i:i + 3] for i in range(0, len(seq) - 2, 3) if "N" not in seq[i:i + 3]
    ))


def frequencia_aminoacidos(proteina: str) -> dict[str, int]:
    """Conta os aminoácidos de uma proteína (ignora '*' e 'X')."""
    return dict(Counter(aa for aa in proteina if aa not in "*X"))


# ---------------------------------------------------------------------------
# Primers
# ---------------------------------------------------------------------------

def temperatura_melting(primer: str) -> float:
    """Estimativa da temperatura de melting (°C) de um oligonucleotídeo.

    - Até 13 nt: regra de Wallace, Tm = 2(A+T) + 4(G+C)
    - Acima de 13 nt: Tm = 64,9 + 41 (G+C - 16,4) / N

    São aproximações clássicas, úteis para triagem de primers.

    >>> temperatura_melting("ATGCATGC")
    24.0
    """
    primer = _exigir_valida(primer)
    at = primer.count("A") + primer.count("T")
    gc = primer.count("G") + primer.count("C")
    n = at + gc
    if n == 0:
        raise ValueError("Primer sem bases A/C/G/T.")
    if n <= 13:
        return float(2 * at + 4 * gc)
    return round(64.9 + 41 * (gc - 16.4) / n, 1)


# ---------------------------------------------------------------------------
# Relatório
# ---------------------------------------------------------------------------

def resumo(seq: str, tamanho_minimo_orf_aa: int = 100) -> dict:
    """Gera um resumo com as principais métricas da sequência."""
    seq = _exigir_valida(seq)
    orfs = encontrar_orfs(seq, tamanho_minimo_orf_aa)
    return {
        "tamanho_pb": len(seq),
        "nucleotideos": contar_nucleotideos(seq),
        "gc_percentual": conteudo_gc(seq),
        "n_orfs": len(orfs),
        "maior_orf_aa": orfs[0].tamanho_aa if orfs else 0,
    }
