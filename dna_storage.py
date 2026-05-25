import struct
from typing import List, Optional, Dict, Tuple
from dna_utils import (bytes_to_bits, bits_to_bytes, bits_to_nucleotide,
                       nucleotide_to_bits, scramble_bits, descramble_bits,
                       find_best_scramble_seed, gc_content, max_homopolymer)
from fountain import FountainEncoder, FountainDecoder, Droplet

class DNAStorageSystem:
    def __init__(self, chunk_size=4, droplet_factor=2.0, scramble_trials=50):
        self.chunk_size = chunk_size
        self.droplet_factor = droplet_factor
        self.scramble_trials = scramble_trials
    def store(self, data, filename="file.bin"):
        encoder = FountainEncoder(data, self.chunk_size, self.droplet_factor)
        droplets = encoder.get_droplets()
        sample = struct.pack('>I', droplets[0].seed) + struct.pack('>H', len(droplets[0].data)) + droplets[0].data + bytes([0])
        bits = bytes_to_bits(sample)
        best_seed = find_best_scramble_seed(bits, self.scramble_trials)
        seqs = []
        for d in droplets:
            payload = struct.pack('>I', d.seed) + struct.pack('>H', len(d.data)) + d.data
            checksum = 0
            for b in payload: checksum ^= b
            payload += bytes([checksum])
            scrambled = scramble_bits(bytes_to_bits(payload), best_seed)
            seqs.append(bits_to_nucleotide(scrambled))
        return {
            'filename': filename,
            'original_size': len(data),
            'num_blocks': encoder.num_blocks,
            'chunk_size': self.chunk_size,
            'scramble_seed': best_seed,
            'num_droplets': len(droplets),
            'dna_sequences': seqs,
            'total_nucleotides': sum(len(s) for s in seqs),
            'gc_content_avg': sum(gc_content(s) for s in seqs)/len(seqs),
            'max_homopolymer_max': max(max_homopolymer(s) for s in seqs)
        }
    def retrieve(self, seqs, meta, error_subst=0, error_del=0, error_loss=0):
        import random
        decoder = FountainDecoder(meta['num_blocks'], meta['chunk_size'])
        stats = {'total_droplets': len(seqs), 'valid':0, 'corrupted':0, 'lost':0}
        for seq in seqs:
            if random.random() < error_loss:
                stats['lost'] += 1
                continue
            bits = nucleotide_to_bits(seq)
            descrambled = descramble_bits(bits, meta['scramble_seed'])
            try:
                payload = bits_to_bytes(descrambled)
            except:
                stats['corrupted'] += 1
                continue
            if len(payload) < 7:
                stats['corrupted'] += 1
                continue
            seed = struct.unpack('>I', payload[0:4])[0]
            data_len = struct.unpack('>H', payload[4:6])[0]
            if len(payload) != 6 + data_len + 1:
                stats['corrupted'] += 1
                continue
            data = payload[6:6+data_len]
            checksum = payload[6+data_len]
            calc = 0
            for b in payload[:6+data_len]: calc ^= b
            if calc != checksum:
                stats['corrupted'] += 1
                continue
            stats['valid'] += 1
            decoder.add_droplet(Droplet(seed, data))
        decoded = decoder.decode()
        if decoded is not None:
            return decoded[:meta['original_size']], stats
        return None, stats