from lazy_loader import LazyLoader 

SEQ_REGION_HEADER = ['seq_region_id', 'name', 'coord_system_id']

loader = LazyLoader("seq_region.txt.gz", SEQ_REGION_HEADER)

