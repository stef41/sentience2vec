#!/usr/bin/env python3
"""Check current data collection status"""

import json
import os
from pathlib import Path

base = Path('/data/users/zacharie/whc/consciousness_study/data/raw')

print('='*70)
print('DATA COLLECTION STATUS')
print('='*70)

# PsychonautWiki
pw_file = base / 'psychonautwiki/substances.json'
if pw_file.exists():
    with open(pw_file) as f:
        pw_data = json.load(f)
    with open(base / 'psychonautwiki/effects.json') as f:
        effects = json.load(f)
    total_effects = sum(len(s.get('effects', [])) for s in pw_data)
    print(f'\n📊 PSYCHONAUTWIKI')
    print(f'   Substances: {len(pw_data)}')
    print(f'   Unique effects: {len(effects)}')
    print(f'   Effect mappings: {total_effects}')

# TripSit
ts_file = base / 'tripsit/key_substances.json'
if ts_file.exists():
    with open(ts_file) as f:
        ts_data = json.load(f)
    combos = sum(len(s.get('combos', {})) for s in ts_data)
    print(f'\n📊 TRIPSIT')
    print(f'   Substances: {len(ts_data)}')
    print(f'   Combo entries: {combos}')

# Erowid
er_dir = base / 'erowid/reports'
if er_dir.exists():
    reports = [f for f in os.listdir(er_dir) if f.endswith('.json')]
    total_chars = 0
    for f in reports:
        with open(er_dir / f) as fp:
            total_chars += len(json.load(fp).get('text', ''))
    print(f'\n📊 EROWID')
    print(f'   Full reports: {len(reports)}')
    print(f'   Text corpus: {total_chars:,} chars ({total_chars/1000000:.1f}MB)')

print('\n' + '='*70)
