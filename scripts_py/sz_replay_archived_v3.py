#!/usr/bin/env python3
"""Extends archived replay v2 with the historical HR/NR Basal 4 Hz cases.
The original v2 engine and all run configurations are preserved.
"""
import sz_replay_archived_v2 as engine

engine.CASES.update({
    'HR_basal4hz': '20250628122534_basal_activity_15_0_glutamate_all',
    'NR_basal4hz': '20251004180523_basal_activity_2_0_glutamate_all',
})

if __name__ == '__main__':
    engine.main()
