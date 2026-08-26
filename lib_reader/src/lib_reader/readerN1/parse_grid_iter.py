from pathlib import Path

import numpy as np

from .parse_grid_data import parse_grid_data_new


def parse_grid_data_iter(file_name, xml_file=None, data_tag='grid1x_wf_packet',
        multi_evt=None, multi_step=None, endian='MSB', crc_check=True,
        skip_et=None, packet_len=1080, chunk_events=50000):
    if skip_et is None:
        skip_et = []

    packet_len = int(packet_len)
    chunk_events = int(chunk_events)
    if packet_len <= 0:
        raise ValueError('packet_len must be positive')
    if chunk_events <= 0:
        raise ValueError('chunk_events must be positive')

    file_path = Path(file_name)
    file_size = file_path.stat().st_size
    if file_size % packet_len != 0:
        raise ValueError(
            f'{file_name} size {file_size} is not aligned to packet_len={packet_len}'
        )

    total_events = file_size // packet_len
    with file_path.open('rb') as fin:
        for chunk_index, event_start in enumerate(range(0, total_events, chunk_events)):
            event_count = min(chunk_events, total_events - event_start)
            raw = fin.read(event_count * packet_len)
            if len(raw) != event_count * packet_len:
                raise IOError(
                    f'failed to read {event_count * packet_len} bytes from {file_name}'
                )

            chunk_data = np.frombuffer(raw, dtype=np.uint8)
            parsed, index = parse_grid_data_new(
                file_name=file_name,
                xml_file=xml_file,
                data_tag=data_tag,
                multi_evt=multi_evt,
                multi_step=multi_step,
                endian=endian,
                crc_check=crc_check,
                data=chunk_data,
                skip_et=skip_et,
                packet_len=packet_len,
            )
            yield {
                'data': parsed,
                'index': index,
                'chunk_index': chunk_index,
                'event_start': event_start,
                'event_count': event_count,
                'total_events': total_events,
            }