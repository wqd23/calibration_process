from lib_reader import single_read10
from .. import util_lib as util
from pathlib import Path
path = Path('data/10B/raw_data/tb_data')

for file in path.iterdir():
    if 'observe' not in file.name:
        continue
    sciExtracted, telExtracted = single_read10(str(file))
    util.pickle_save({'sci':sciExtracted, 'tel': telExtracted},f'data/10B/raw_data/tb_cache/10B_{file.stem}.pkl')