# -*- coding:utf-8 -*-
"""Per-version explicit workflows.

Each version module encodes the historical selection rules and any genuinely
version-specific orchestration, so the full flow can be read from that single
file.  Version selection happens once here (in ``registry``): after a version
module is imported no ``if version == ...`` branch should appear.
"""
