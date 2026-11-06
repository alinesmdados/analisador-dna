"""Testes automatizados do analisador_dna (execute com: pytest)."""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import analisador_dna as ad  # noqa: E402


def test_ler_fasta_multiplas_sequencias():
    texto = ">seq1 primeira\nATG\nCCC\n\n>seq2\nggg tt\n"
    regs = ad.ler_fasta_texto(texto)
    assert [r.identificador for r in regs] == ["seq1", "seq2"]
    assert regs[0].descricao == "primeira"
    assert regs[0].sequencia == "ATGCCC"
    assert regs[1].sequencia == "GGGTT"


def test_fasta_sem_cabecalho_gera_erro():
    with pytest.raises(ValueError):
        ad.ler_fasta_texto("ATGC\n")


def test_validacao():
    assert ad.validar_sequencia("ACGTN")
    assert not ad.validar_sequencia("ACGU")
    assert not ad.validar_sequencia("")
    with pytest.raises(ValueError):
        ad.conteudo_gc("ATGX")


def test_contagem_e_gc():
    assert ad.contar_nucleotideos("AACGTT") == {"A": 2, "C": 1, "G": 1, "T": 2}
    assert ad.conteudo_gc("GGCC") == 100.0
    assert ad.conteudo_gc("ATAT") == 0.0
    assert ad.conteudo_gc("GCNN") == 100.0  # N fora do denominador


def test_gc_janela_deslizante():
    seq = "A" * 50 + "G" * 50
    pontos = ad.gc_janela_deslizante(seq, janela=50, passo=50)
    assert pontos == [(25, 0.0), (75, 100.0)]


def test_complemento_reverso_e_involutivo():
    seq = "ATGCGTACGTTAGC"
    assert ad.complemento_reverso("ATGCCG") == "CGGCAT"
    assert ad.complemento_reverso(ad.complemento_reverso(seq)) == seq


def test_transcricao_e_traducao():
    assert ad.transcrever("ATGCTT") == "AUGCUU"
    assert ad.traduzir("ATGGCCTAA") == "MA*"
    assert ad.traduzir("ATGGCCTAAGGG", parar_no_stop=True) == "MA"
    assert ad.traduzir("ATGNNN") == "MX"


def test_codigo_genetico_completo():
    assert len(ad.CODIGO_GENETICO) == 64
    assert ad.CODIGO_GENETICO["ATG"] == "M"
    assert ad.CODIGO_GENETICO["TGG"] == "W"
    assert {c for c, aa in ad.CODIGO_GENETICO.items() if aa == "*"} == ad.CODONS_PARADA


def test_orf_fita_direta():
    # 5 nt de enchimento + ATG + 10 códons GCT + TAA
    seq = "CCCCC" + "ATG" + "GCT" * 10 + "TAA" + "CCCCC"
    orfs = ad.encontrar_orfs(seq, tamanho_minimo_aa=5)
    direta = [o for o in orfs if o.fita == "+"]
    assert len(direta) == 1
    o = direta[0]
    assert o.proteina == "M" + "A" * 10
    assert (o.inicio, o.fim) == (6, 41)
    assert seq[o.inicio - 1:o.fim].startswith("ATG")
    assert seq[o.inicio - 1:o.fim].endswith("TAA")


def test_orf_fita_reversa_coordenadas():
    gene = "ATG" + "GCT" * 10 + "TAA"
    seq = "CCCCC" + ad.complemento_reverso(gene) + "CCCCC"
    orfs = [o for o in ad.encontrar_orfs(seq, tamanho_minimo_aa=5) if o.fita == "-"]
    assert len(orfs) == 1
    o = orfs[0]
    trecho = seq[o.inicio - 1:o.fim]
    assert ad.complemento_reverso(trecho) == gene


def test_motivo_palindromico_nao_duplica():
    # GAATTC (sítio da EcoRI) é palíndromo: não deve aparecer duas vezes
    assert ad.encontrar_motivo("AAGAATTCAA", "GAATTC") == [(3, "+")]


def test_motivo_fita_reversa():
    seq = "AAAACCGT"  # ACGG na fita reversa = CCGT na direta
    assert ad.encontrar_motivo(seq, "ACGG") == [(5, "-")]


def test_frequencia_codons():
    assert ad.frequencia_codons("ATGATGTAA") == {"ATG": 2, "TAA": 1}


def test_temperatura_melting():
    assert ad.temperatura_melting("ATGCATGC") == 24.0
    tm = ad.temperatura_melting("AGCGGATAACAATTTCACACAGGA")  # primer M13 reverse
    assert 50 < tm < 60


def test_resumo():
    r = ad.resumo("ATG" + "GCT" * 120 + "TAA", tamanho_minimo_orf_aa=100)
    assert r["tamanho_pb"] == 366
    assert r["n_orfs"] >= 1
    assert r["maior_orf_aa"] == 121
