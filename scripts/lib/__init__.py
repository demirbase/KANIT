#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared library for the AMR Prediction pipeline (SCALE_MLOPS_PLAN.md §5).

Code shared by the numbered scripts.

Submodules:
    registry      — organisms.yaml / antibiotics.yaml access (single source)
    config        — global config loader + {organism}/{antibiotic} path resolver
    chunking      — get_y_chunk (contiguous label slicing)
    io_utils      — run_command (shlex-based, never shell=True)
    run_metadata  — git hash / version capture, run_id generation (MLOps)
"""
