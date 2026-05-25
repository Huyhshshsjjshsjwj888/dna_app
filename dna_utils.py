import random, hashlib, math
from collections import Counter
from typing import Dict, Tuple, List

MAPPING = {'00':'A','01':'T','10':'C','11':'G'}
REVERSE_MAPPING = {'A':'00','T':'01','C':'10','G':'11'}

# --- CƠ BẢN ---
def bytes_to_bits(data): return ''.join(f'{b:08b}' for b in data) if data else ""
def bits_to_bytes(bits):
    if not bits: return b""
    bits += '0' * ((8 - len(bits)%8)%8)
    return bytes(int(bits[i:i+8],2) for i in range(0,len(bits),8))
def bits_to_nucleotide(bits, mapping=None):
    if mapping is None: mapping = MAPPING
    if not bits: return ""
    if len(bits)%2: bits += '0'
    return ''.join(mapping[bits[i:i+2]] for i in range(0,len(bits),2))
def nucleotide_to_bits(seq, rev=None):
    if rev is None: rev = REVERSE_MAPPING
    return ''.join(rev[n] for n in seq) if seq else ""
def scramble_bits(bits, seed):
    if not bits: return ""
    r = random.Random(seed)
    return ''.join(str(int(b)^r.randint(0,1)) for b in bits)
def descramble_bits(bits, seed): return scramble_bits(bits, seed)
def gc_content(seq): return (seq.count('C')+seq.count('G'))/len(seq) if seq else 0.0
def max_homopolymer(seq):
    if not seq: return 0
    m=c=1
    for i in range(1,len(seq)):
        if seq[i]==seq[i-1]: c+=1; m=max(m,c)
        else: c=1
    return m
def find_best_scramble_seed(bits, n=50):
    if not bits: return 0
    best=0; best_score=float('inf')
    for s in range(n):
        seq=bits_to_nucleotide(scramble_bits(bits,s))
        score=abs(gc_content(seq)-0.5)*100+max(0,max_homopolymer(seq)-3)*10
        if score<best_score: best_score=score; best=s
    return best
def simulate_errors(seq, sub=0, dl=0, ins=0):
    if not seq: return ""
    nuc=['A','T','C','G']; r=[]
    for ch in seq:
        if random.random()<dl: continue
        if random.random()<sub: ch=random.choice([x for x in nuc if x!=ch])
        r.append(ch)
        if random.random()<ins: r.append(random.choice(nuc))
    return ''.join(r)

# --- HUFFMAN ---
class HuffmanNode:
    def __init__(self, char=None, freq=0): self.char=char; self.freq=freq; self.left=None; self.right=None
def build_huffman_tree(data):
    if not data: return None
    nodes=[HuffmanNode(b,c) for b,c in Counter(data).items()]
    while len(nodes)>1:
        nodes.sort(key=lambda x:x.freq)
        l=nodes.pop(0); r=nodes.pop(0)
        p=HuffmanNode(freq=l.freq+r.freq); p.left=l; p.right=r
        nodes.append(p)
    return nodes[0]
def generate_huffman_codes(node, code="", codes=None):
    if codes is None: codes={}
    if node is None: return codes
    if node.char is not None: codes[node.char]=code or "0"; return codes
    generate_huffman_codes(node.left,code+"0",codes)
    generate_huffman_codes(node.right,code+"1",codes)
    return codes
def huffman_compress(data):
    if not data: return "",{},0
    codes=generate_huffman_codes(build_huffman_tree(data))
    return ''.join(codes[b] for b in data), codes, len(data)
def huffman_decompress(bits, codes, orig_len):
    if not bits or not codes: return b""
    rev={v:k for k,v in codes.items()}
    res=bytearray(); cur=""
    for b in bits:
        cur+=b
        if cur in rev: res.append(rev[cur]); cur=""; 
        if len(res)>=orig_len: break
    return bytes(res)

# --- AES ---
def aes_encrypt(data, pwd):
    if not pwd: return data
    key=hashlib.sha256(pwd.encode()).digest()
    return bytes(b^key[i%len(key)] for i,b in enumerate(data))
def aes_decrypt(data, pwd): return aes_encrypt(data, pwd)

# --- GENETIC ALGORITHM ---
class GeneticOptimizer:
    def __init__(self, bits, pop=20, gen=30): self.bits=bits; self.pop=pop; self.gen=gen
    def fitness(self, seed):
        seq=bits_to_nucleotide(scramble_bits(self.bits, seed))
        return -(abs(gc_content(seq)-0.5)*100+max(0,max_homopolymer(seq)-3)*15)
    def crossover(self, p1, p2):
        b1=format(p1,'032b'); b2=format(p2,'032b')
        pt=random.randint(1,30)
        return int(b1[:pt]+b2[pt:],2)%1000000
    def mutate(self, s, rate=0.1): return random.randint(0,999999) if random.random()<rate else s
    def optimize(self):
        pop=[random.randint(0,999999) for _ in range(self.pop)]
        for _ in range(self.gen):
            scores=[(s,self.fitness(s)) for s in pop]
            scores.sort(key=lambda x:x[1], reverse=True)
            top=[s for s,_ in scores[:self.pop//2]]
            new=top.copy()
            while len(new)<self.pop:
                p1,p2=random.choice(top),random.choice(top)
                new.append(self.mutate(self.crossover(p1,p2)))
            pop=new
        return max(pop, key=lambda s:self.fitness(s))

# --- 3D ---
def generate_dna_3d_coordinates(seq):
    coords=[]; r=1.0; h=0.34; a=36
    comp={'A':'T','T':'A','C':'G','G':'C'}
    for i,b in enumerate(seq[:100]):
        ang=math.radians(i*a); z=i*h
        coords.append({'x':round(r*math.cos(ang),3),'y':round(r*math.sin(ang),3),'z':round(z,3),'nucleotide':b,'strand':1,'position':i})
        coords.append({'x':round(r*math.cos(ang+math.pi),3),'y':round(r*math.sin(ang+math.pi),3),'z':round(z,3),'nucleotide':comp[b],'strand':2,'position':i})
    return coords

# --- PCR ---
def simulate_pcr(seq, cycles=3, err=0.001):
    copies=[seq]; total_err=0
    for _ in range(cycles):
        new=[]
        for c in copies:
            m=simulate_errors(c, err, err/10, err/10)
            new.append(m)
            if m!=c: total_err+=1
        copies.extend(new)
    return copies, total_err

# --- ENZYME & PRESETS ---
ENZYME_INFO = {
    "Taq":{"base_subst":0.01,"base_del":0.001,"base_ins":0.001,"temp_factor":0.1,"description":"Taq"},
    "Phusion":{"base_subst":0.0005,"base_del":0.0002,"base_ins":0.0002,"temp_factor":0.05,"description":"Phusion"},
    "Q5":{"base_subst":0.0001,"base_del":0.00005,"base_ins":0.00005,"temp_factor":0.02,"description":"Q5"}
}
ERROR_PRESETS = {
    "Tổng hợp chính xác cao":{"subst":0.001,"del":0.0005,"ins":0.0005},
    "PCR Taq":{"subst":0.01,"del":0.001,"ins":0.001},
    "Nanopore":{"subst":0.05,"del":0.03,"ins":0.02},
    "Illumina":{"subst":0.001,"del":0.0001,"ins":0.0001},
    "Suy thoái tự nhiên":{"subst":0.02,"del":0.005,"ins":0.002}
}
def get_error_rates_from_enzyme(name, temp=55.0):
    e=ENZYME_INFO.get(name, ENZYME_INFO["Taq"])
    delta=temp-55.0; factor=1.0+e["temp_factor"]*delta
    factor=max(0.5, min(2.0, factor))
    return {"subst":e["base_subst"]*factor, "del":e["base_del"]*factor, "ins":e["base_ins"]*factor}
def simulate_errors_advanced(seq, mode="preset", preset_name=None, enzyme=None, temperature=55.0):
    if mode=="preset" and preset_name:
        r=ERROR_PRESETS.get(preset_name)
        if r: return simulate_errors(seq, r["subst"], r["del"], r["ins"])
    elif mode=="enzyme" and enzyme:
        r=get_error_rates_from_enzyme(enzyme, temperature)
        return simulate_errors(seq, r["subst"], r["del"], r["ins"])
    return seq