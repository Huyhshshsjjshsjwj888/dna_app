import random, struct
from dataclasses import dataclass
from typing import List, Optional

@dataclass
class Droplet:
    seed: int
    data: bytes

class FountainEncoder:
    def __init__(self, data: bytes, chunk_size=4, factor=2.0, seed_base=12345):
        self.data = data
        self.chunk_size = chunk_size
        pad = (chunk_size - len(data) % chunk_size) % chunk_size
        self.padded = data + b'\x00'*pad
        self.orig_size = len(data)
        self.blocks = [self.padded[i:i+chunk_size] for i in range(0, len(self.padded), chunk_size)]
        self.num_blocks = len(self.blocks)
        self.num_droplets = max(int(self.num_blocks*factor), self.num_blocks+5)
        self.seed_base = seed_base
        self.droplets = []
        self._gen()
    def _gen(self):
        for i in range(self.num_droplets):
            seed = self.seed_base + i
            rng = random.Random(seed)
            idxs = [idx for idx in range(self.num_blocks) if rng.random()<0.5]
            if not idxs: idxs = [rng.randint(0,self.num_blocks-1)]
            xor = bytearray(self.chunk_size)
            for idx in idxs:
                block = self.blocks[idx]
                for j in range(self.chunk_size): xor[j] ^= block[j]
            self.droplets.append(Droplet(seed, bytes(xor)))
    def get_droplets(self): return self.droplets

class FountainDecoder:
    def __init__(self, num_blocks, chunk_size, seed_base=12345):
        self.num_blocks = num_blocks
        self.chunk_size = chunk_size
        self.seed_base = seed_base
        self.droplets = []
        self.solved = [None]*num_blocks
        self.solved_cnt = 0
    def add_droplet(self, d): self.droplets.append(d)
    def decode(self):
        infos = []
        for d in self.droplets:
            rng = random.Random(d.seed)
            idxs = set()
            for idx in range(self.num_blocks):
                if rng.random()<0.5: idxs.add(idx)
            if not idxs: idxs.add(rng.randint(0,self.num_blocks-1))
            infos.append({'seed':d.seed, 'data':bytearray(d.data), 'indices':idxs})
        changed = True
        while changed and self.solved_cnt < self.num_blocks:
            changed = False
            for info in infos:
                unresolved = [idx for idx in info['indices'] if self.solved[idx] is None]
                if len(unresolved)==1:
                    target = unresolved[0]
                    mask = bytearray(self.chunk_size)
                    for idx in info['indices']:
                        if idx!=target and self.solved[idx] is not None:
                            block = self.solved[idx]
                            for j in range(self.chunk_size): mask[j] ^= block[j]
                    solved_block = bytearray(self.chunk_size)
                    for j in range(self.chunk_size): solved_block[j] = info['data'][j] ^ mask[j]
                    self.solved[target] = bytes(solved_block)
                    self.solved_cnt += 1
                    changed = True
                    for o in infos:
                        if target in o['indices']:
                            for j in range(self.chunk_size): o['data'][j] ^= solved_block[j]
                            o['indices'].discard(target)
                    break
        return b''.join(self.solved) if self.solved_cnt==self.num_blocks else None