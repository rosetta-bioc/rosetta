"""Tests for the opt-in R code-generation log."""

from rosetta import codegen


def setup_function():
    codegen.disable()
    codegen.clear()


def teardown_function():
    codegen.disable()
    codegen.clear()


def test_emit_is_ignored_when_codegen_is_disabled(capsys):
    codegen._emit("dds <- DESeq(dds)")

    assert codegen.last() == ""
    assert capsys.readouterr().out == ""


def test_emit_records_and_prints_when_codegen_is_enabled(capsys):
    codegen.enable()

    codegen._emit("dds <- DESeq(dds)")

    assert codegen.last() == "dds <- DESeq(dds)"
    assert "dds <- DESeq(dds)" in capsys.readouterr().out
