from flask import Flask, request, jsonify, send_file
from dna_utils import (simulate_errors, get_error_rates_from_enzyme,
                       huffman_compress, huffman_decompress,
                       aes_encrypt, aes_decrypt,
                       GeneticOptimizer,
                       generate_dna_3d_coordinates, simulate_pcr,
                       ENZYME_INFO, ERROR_PRESETS,
                       bits_to_bytes, bytes_to_bits)
from dna_storage import DNAStorageSystem
import io, random, json
from datetime import datetime

app = Flask(__name__)
stored = {}
history = []

HTML = """<!DOCTYPE html>
<html lang="vi">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0, maximum-scale=1.0">
<title>🧬 DNA Storage Pro</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three.js@0.160.0/build/three.min.js"></script>
<style>
:root{--bg:#0a0e1a;--surface:#111827;--surface2:#1a2235;--accent:#00ff9d;--accent2:#00bfff;--danger:#ff4d6d;--warn:#ffd166;--text:#e2e8f0;--muted:#64748b;--border:#1e293b}
*{box-sizing:border-box;margin:0;padding:0}
body{background:var(--bg);color:var(--text);font-family:'Syne',sans-serif;min-height:100vh;padding:1rem}
.container{max-width:720px;margin:0 auto}
header{text-align:center;padding:1.5rem 0}
header h1{font-size:2rem;font-weight:800;background:linear-gradient(135deg,var(--accent),var(--accent2));-webkit-background-clip:text;-webkit-text-fill-color:transparent;background-clip:text}
.tabs{display:flex;gap:4px;overflow-x:auto;padding-bottom:.5rem;scrollbar-width:none}
.tabs::-webkit-scrollbar{display:none}
.tab-btn{padding:.6rem .8rem;background:var(--surface);border:1px solid var(--border);border-radius:8px;color:var(--muted);font-size:.8rem;font-weight:600;cursor:pointer;white-space:nowrap;transition:.2s}
.tab-btn.active{background:var(--surface2);border-color:var(--accent);color:var(--accent)}
.panel{display:none;padding:1rem 0}
.panel.active{display:block}
.card{background:var(--surface);border:1px solid var(--border);border-radius:16px;padding:1.25rem;margin-bottom:1rem}
.card-title{font-size:1.1rem;font-weight:800;margin-bottom:1rem}
textarea,input[type=file],select,input[type=password]{width:100%;background:var(--surface2);border:1px solid var(--border);border-radius:10px;color:var(--text);font-family:'Space Mono',monospace;font-size:.9rem;padding:.75rem;resize:vertical;outline:none;margin-bottom:1rem}
textarea{min-height:120px}
.btn{width:100%;padding:1rem;border:none;border-radius:12px;font-family:'Syne',sans-serif;font-size:1rem;font-weight:700;cursor:pointer;transition:.2s;margin-top:.5rem;letter-spacing:.5px}
.btn-primary{background:linear-gradient(135deg,var(--accent),#00cc7a);color:#0a0e1a}
.btn-secondary{background:linear-gradient(135deg,var(--accent2),#0080cc);color:#0a0e1a}
.btn-danger{background:linear-gradient(135deg,var(--danger),#cc0033);color:white}
.btn-warning{background:linear-gradient(135deg,var(--warn),#ff9900);color:#0a0e1a}
.btn:disabled{opacity:.4;cursor:not-allowed}
.spinner{display:none;text-align:center;padding:1rem;color:var(--accent);font-family:'Space Mono',monospace}
.spinner.show{display:block}
.spinner::before{content:"";display:block;width:30px;height:30px;border:3px solid var(--border);border-top-color:var(--accent);border-radius:50%;margin:0 auto .5rem;animation:spin .8s linear infinite}
@keyframes spin{to{transform:rotate(360deg)}}
.alert{padding:.8rem;border-radius:10px;margin-top:.5rem;font-family:'Space Mono',monospace;display:none}
.alert.show{display:block}
.alert-success{background:rgba(0,255,157,.1);border:1px solid var(--accent);color:var(--accent)}
.alert-error{background:rgba(255,77,109,.1);border:1px solid var(--danger);color:var(--danger)}
.alert-info{background:rgba(0,191,255,.1);border:1px solid var(--accent2);color:var(--accent2)}
.alert-warning{background:rgba(255,209,102,.1);border:1px solid var(--warn);color:var(--warn)}
.metrics{display:grid;grid-template-columns:1fr 1fr;gap:.5rem;margin-top:1rem}
.metric{background:var(--surface2);border:1px solid var(--border);border-radius:10px;padding:.7rem;text-align:center}
.metric-value{font-size:1.2rem;font-weight:800;color:var(--accent);font-family:'Space Mono',monospace}
.metric-label{font-size:.7rem;color:var(--muted)}
.slider-row{display:flex;justify-content:space-between;margin-bottom:.3rem}
input[type=range]{width:100%;accent-color:var(--accent)}
canvas{background:var(--surface2);border-radius:10px;margin-top:1rem}
.min-droplet-info,.max-loss-info{background:rgba(0,191,255,.08);border:1px solid var(--accent2);border-radius:8px;padding:.8rem;margin:1rem 0;font-family:'Space Mono',monospace;font-size:.8rem;color:var(--accent2)}
#dna3DContainer{width:100%;height:300px;background:var(--surface2);border-radius:10px;margin-top:1rem}
.history-table{width:100%;border-collapse:collapse;margin-top:1rem;font-size:.8rem}
.history-table th,.history-table td{padding:.5rem;border:1px solid var(--border);text-align:center}
.history-table th{background:var(--surface2);color:var(--accent)}
.radio-group{display:flex;gap:.5rem;margin-bottom:1rem}
.radio-option{flex:1;position:relative}
.radio-option input{position:absolute;opacity:0}
.radio-option label{display:block;padding:.7rem;background:var(--surface2);border:2px solid var(--border);border-radius:10px;text-align:center;cursor:pointer;font-weight:600;font-size:.9rem;transition:.2s}
.radio-option input:checked+label{border-color:var(--accent);color:var(--accent);background:rgba(0,255,157,.08)}
.file-section{display:none}
</style>
</head>
<body>
<div class="container">
<header><h1>🧬 DNA Storage Pro</h1><p style="color:var(--muted);font-size:.75rem">Mô phỏng lưu trữ dữ liệu bằng DNA - Phiên bản Quốc gia</p></header>
<div class="tabs" id="tabMenu">
<button class="tab-btn active" onclick="switchTab('encode',this)">📂 Mã hoá</button>
<button class="tab-btn" onclick="switchTab('decode',this)">🔍 Giải mã</button>
<button class="tab-btn" onclick="switchTab('error',this)">⚠️ Lỗi</button>
<button class="tab-btn" onclick="switchTab('compare',this)">🔬 So sánh</button>
<button class="tab-btn" onclick="switchTab('reality',this)">🧪 Thực tế ảo</button>
<button class="tab-btn" onclick="switchTab('pcr',this)">🧬 PCR ảo</button>
<button class="tab-btn" onclick="switchTab('3d',this)">🧬 3D DNA</button>
<button class="tab-btn" onclick="switchTab('history',this)">📜 Lịch sử</button>
<button class="tab-btn" onclick="switchTab('life',this)">📊 Tuổi thọ</button>
</div>

<!-- TAB MÃ HOÁ -->
<div id="tab-encode" class="panel active">
<div class="card">
<div class="card-title">1. Nhập dữ liệu</div>
<div class="radio-group">
<div class="radio-option"><input type="radio" name="inputMode" id="modeText" value="text" checked onchange="toggleMode('text')"><label for="modeText">✏️ Văn bản</label></div>
<div class="radio-option"><input type="radio" name="inputMode" id="modeFile" value="file" onchange="toggleMode('file')"><label for="modeFile">📁 Tải file</label></div>
</div>
<div id="textSection">
<textarea id="textInput" placeholder="Gõ nội dung..."></textarea>
</div>
<div class="file-section" id="fileSection">
<input type="file" id="fileInput" onchange="onFileSelected()">
<div class="alert alert-info" id="fileInfo"></div>
</div>
<label>🔐 Mật khẩu AES (tùy chọn):</label>
<input type="password" id="aesPassword" placeholder="Để trống nếu không mã hóa">
<div class="radio-group">
<div class="radio-option"><input type="radio" name="optMode" id="optAI" value="ai" checked><label for="optAI">🤖 AI Tối ưu</label></div>
<div class="radio-option"><input type="radio" name="optMode" id="optGA" value="ga"><label for="optGA">🧬 Genetic Algorithm</label></div>
</div>
<label>🧬 Nén Huffman:</label>
<div class="radio-group">
<div class="radio-option"><input type="radio" name="huffmanMode" id="huffmanOn" value="on" checked><label for="huffmanOn">✅ Bật</label></div>
<div class="radio-option"><input type="radio" name="huffmanMode" id="huffmanOff" value="off"><label for="huffmanOff">❌ Tắt</label></div>
</div>
<button class="btn btn-primary" id="encodeBtn" onclick="encodeData()">🧬 Mã hoá DNA</button>
<div class="spinner" id="encodeSpinner">Đang mã hoá...</div>
<div class="alert" id="encodeAlert"></div>
</div>
<div class="card" id="resultCard" style="display:none">
<div class="card-title">🔬 Kết quả</div>
<div class="metrics" id="metricsGrid"></div>
<div class="alert alert-info" id="dnaPreview" style="font-size:.75rem;word-break:break-all"></div>
<div class="min-droplet-info" id="minDropletInfo" style="display:none"></div>
<div class="max-loss-info" id="maxLossInfo" style="display:none"></div>
</div>
</div>

<!-- TAB GIẢI MÃ -->
<div id="tab-decode" class="panel">
<div class="card">
<div class="card-title">2. Giải mã & Phục hồi</div>
<div id="decodeStatus" class="alert alert-info show">⚠️ Chưa có dữ liệu</div>
<label>🔐 Mật khẩu (nếu có):</label>
<input type="password" id="decodePassword" placeholder="Nhập mật khẩu">
<button class="btn btn-secondary" id="decodeBtn" onclick="decodeData()" disabled>🔬 Giải mã</button>
<div class="spinner" id="decodeSpinner"></div>
<div class="alert" id="decodeAlert"></div>
<textarea id="recoveredText" readonly style="min-height:100px;color:var(--accent);display:none"></textarea>
<button class="btn btn-secondary" onclick="downloadFile()" style="display:none" id="downloadBtn">⬇️ Tải file</button>
</div>
</div>

<!-- TAB LỖI -->
<div id="tab-error" class="panel">
<div class="card">
<div class="card-title">⚠️ Mô phỏng lỗi</div>
<div id="errorStatus" class="alert alert-info show">⚠️ Hãy mã hoá trước</div>
<div id="minDropletInfoError" class="min-droplet-info" style="display:none"></div>
<div id="maxLossInfoError" class="max-loss-info" style="display:none"></div>
<div class="slider-row"><span>Đột biến thay thế</span><span id="substVal">1%</span></div>
<input type="range" min="0" max="50" value="1" oninput="document.getElementById('substVal').textContent=this.value+'%'" id="substRate">
<div class="slider-row"><span>Xoá nucleotide</span><span id="delVal">1%</span></div>
<input type="range" min="0" max="30" value="1" oninput="document.getElementById('delVal').textContent=this.value+'%'" id="delRate">
<div class="slider-row"><span>Mất droplet</span><span id="lossVal">5%</span></div>
<input type="range" min="0" max="50" value="5" oninput="document.getElementById('lossVal').textContent=this.value+'%'" id="lossRate">
<button class="btn btn-danger" id="errorBtn" onclick="simulateError()" disabled>💥 Phá huỷ & Phục hồi</button>
<div class="spinner" id="errorSpinner"></div>
<div class="alert" id="errorAlert"></div>
<div class="metrics" id="errorMetrics" style="display:none"></div>
</div>
</div>

<!-- TAB SO SÁNH -->
<div id="tab-compare" class="panel">
<div class="card">
<div class="card-title">🔬 So sánh hiệu năng sửa lỗi</div>
<div id="compareStatus" class="alert alert-info show">⚠️ Hãy mã hoá trước</div>
<div class="slider-row"><span>Tỉ lệ mất droplet</span><span id="compareLossVal">10%</span></div>
<input type="range" min="0" max="80" value="10" oninput="document.getElementById('compareLossVal').textContent=this.value+'%'" id="compareLossRate">
<button class="btn btn-primary" id="compareBtn" onclick="runComparison()" disabled>📊 So sánh</button>
<canvas id="compareChart" height="250"></canvas>
</div>
</div>

<!-- TAB THỰC TẾ ẢO -->
<div id="tab-reality" class="panel">
<div class="card">
<div class="card-title">🧪 Mô phỏng enzyme & nhiệt độ</div>
<div id="realityStatus" class="alert alert-info show">⚠️ Hãy mã hoá trước</div>
<label>Enzyme:</label>
<select id="enzymeSelect">
<option value="Taq">Taq polymerase</option>
<option value="Phusion" selected>Phusion</option>
<option value="Q5">Q5</option>
</select>
<label>Nhiệt độ annealing (°C):</label>
<div class="slider-row"><span>Nhiệt độ</span><span id="tempVal">55°C</span></div>
<input type="range" min="45" max="65" value="55" oninput="document.getElementById('tempVal').textContent=this.value+'°C'" id="tempSlider">
<button class="btn btn-danger" id="realityBtn" onclick="simulateReality()" disabled>🧪 Mô phỏng lỗi enzyme</button>
<div class="spinner" id="realitySpinner"></div>
<div class="alert" id="realityAlert"></div>
<div class="metrics" id="realityMetrics" style="display:none"></div>
</div>
</div>

<!-- TAB PCR ẢO -->
<div id="tab-pcr" class="panel">
<div class="card">
<div class="card-title">🧬 PCR ảo</div>
<div id="pcrStatus" class="alert alert-info show">⚠️ Hãy mã hoá trước</div>
<label>Số chu kỳ:</label>
<div class="slider-row"><span>Chu kỳ</span><span id="pcrCyclesVal">3</span></div>
<input type="range" min="1" max="10" value="3" oninput="document.getElementById('pcrCyclesVal').textContent=this.value" id="pcrCycles">
<button class="btn btn-warning" id="pcrBtn" onclick="runPCR()" disabled>🧬 Chạy PCR</button>
<div class="spinner" id="pcrSpinner"></div>
<div class="alert" id="pcrAlert"></div>
<canvas id="pcrChart" height="200"></canvas>
</div>
</div>

<!-- TAB 3D DNA -->
<div id="tab-3d" class="panel">
<div class="card">
<div class="card-title">🧬 Mô hình 3D DNA</div>
<div id="dna3dStatus" class="alert alert-info show">⚠️ Hãy mã hoá trước để xem cấu trúc</div>
<div id="dna3dContainer"></div>
<p style="font-size:.75rem;color:var(--muted)">🖱️ Kéo để xoay</p>
</div>
</div>

<!-- TAB LỊCH SỬ -->
<div id="tab-history" class="panel">
<div class="card">
<div class="card-title">📜 Lịch sử mã hoá</div>
<div id="historyContent"><p style="color:var(--muted)">Chưa có lịch sử</p></div>
<button class="btn btn-primary" onclick="loadHistory()">🔄 Tải lại</button>
</div>
</div>

<!-- TAB TUỔI THỌ -->
<div id="tab-life" class="panel">
<div class="card">
<div class="card-title">📊 Tuổi thọ lưu trữ</div>
<div id="lifeChart"></div>
</div>
</div>
</div>

<script>
let selectedFile=null, chartInstance=null, pcrChartInstance=null;
let threeScene, threeRenderer, threeCamera;
function switchTab(n,b){document.querySelectorAll('.panel').forEach(p=>p.classList.remove('active'));document.querySelectorAll('.tab-btn').forEach(x=>x.classList.remove('active'));document.getElementById('tab-'+n).classList.add('active');b.classList.add('active')}
function showAlert(id,msg,type){const e=document.getElementById(id);e.textContent=msg;e.className='alert alert-'+type+' show'}
function setLoading(id,show){document.getElementById(id).className='spinner'+(show?' show':'')}
function toggleMode(m){document.getElementById('textSection').style.display=m==='text'?'block':'none';document.getElementById('fileSection').style.display=m==='file'?'block':'none'}
function onFileSelected(){const f=document.getElementById('fileInput').files[0];if(f){selectedFile=f;document.getElementById('fileInfo').textContent='✅ '+f.name+' ('+f.size.toLocaleString()+' bytes)';document.getElementById('fileInfo').className='alert alert-success show'}}

async function encodeData(){
const mode=document.querySelector('input[name="inputMode"]:checked').value;
const optMode=document.querySelector('input[name="optMode"]:checked').value;
const huffmanMode=document.querySelector('input[name="huffmanMode"]:checked').value;
const password=document.getElementById('aesPassword').value;
setLoading('encodeSpinner',true);
try{
let res;
if(mode==='text'){
const text=document.getElementById('textInput').value.trim();
if(!text){showAlert('encodeAlert','⚠️ Nhập văn bản','error');setLoading('encodeSpinner',false);return}
res=await fetch('/encode_text',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({text,opt_mode:optMode,huffman:huffmanMode==='on',password})});
}else{
if(!selectedFile){showAlert('encodeAlert','⚠️ Chọn file','error');setLoading('encodeSpinner',false);return}
const fd=new FormData();fd.append('file',selectedFile);fd.append('opt_mode',optMode);fd.append('huffman',huffmanMode==='on');fd.append('password',password);
res=await fetch('/encode_file',{method:'POST',body:fd});
}
const data=await res.json();
if(data.error){showAlert('encodeAlert','❌ '+data.error,'error');setLoading('encodeSpinner',false);return}
showResult(data);
}catch(e){showAlert('encodeAlert','❌ Lỗi kết nối','error')}
finally{setLoading('encodeSpinner',false)}
}

function showResult(data){
showAlert('encodeAlert','✅ Mã hoá thành công!','success');
document.getElementById('resultCard').style.display='block';
document.getElementById('metricsGrid').innerHTML=`<div class="metric"><div class="metric-value">${data.num_droplets}</div><div class="metric-label">Droplet</div></div><div class="metric"><div class="metric-value">${data.total_nucleotides.toLocaleString()}</div><div class="metric-label">Nucleotide</div></div><div class="metric"><div class="metric-value">${(data.gc_content_avg*100).toFixed(1)}%</div><div class="metric-label">GC</div></div><div class="metric"><div class="metric-value">${data.max_homopolymer_max}</div><div class="metric-label">Max lặp</div></div>`;
document.getElementById('dnaPreview').textContent='🧬 '+data.dna_preview;
document.getElementById('dnaPreview').className='alert alert-info show';
['decodeBtn','errorBtn','compareBtn','realityBtn','pcrBtn'].forEach(id=>document.getElementById(id).disabled=false);
['decodeStatus','errorStatus','compareStatus','realityStatus','pcrStatus','dna3dStatus'].forEach(id=>{document.getElementById(id).textContent='✅ Sẵn sàng';document.getElementById(id).className='alert alert-success show'});
if(data.min_droplets){
document.getElementById('minDropletInfo').style.display='block';
document.getElementById('minDropletInfo').textContent='🔍 Cần ít nhất '+data.min_droplets+' droplet để phục hồi.';
document.getElementById('maxLossInfo').style.display='block';
const maxLoss=((data.num_droplets-data.min_droplets)/data.num_droplets*100).toFixed(1);
document.getElementById('maxLossInfo').textContent='⚠️ Tỉ lệ mất tối đa: '+maxLoss+'%';
// Cập nhật cho tab lỗi
document.getElementById('minDropletInfoError').style.display='block';
document.getElementById('minDropletInfoError').textContent='🔍 Cần ít nhất '+data.min_droplets+' droplet hợp lệ.';
document.getElementById('maxLossInfoError').style.display='block';
document.getElementById('maxLossInfoError').textContent='⚠️ Tỉ lệ mất tối đa: '+maxLoss+'%';
}
window._dna3dData = data.dna_3d;
}

async function decodeData(){
const password=document.getElementById('decodePassword').value;
setLoading('decodeSpinner',true);
try{
const res=await fetch('/decode',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({password})});
const data=await res.json();
if(data.error){showAlert('decodeAlert','❌ '+data.error,'error');setLoading('decodeSpinner',false);return}
showAlert('decodeAlert','✅ Phục hồi thành công!','success');
if(data.text){document.getElementById('recoveredText').style.display='block';document.getElementById('recoveredText').value=data.text;document.getElementById('downloadBtn').style.display='block'}
}catch(e){showAlert('decodeAlert','❌ Lỗi','error')}
finally{setLoading('decodeSpinner',false)}
}
function downloadFile(){window.location.href='/download'}

async function simulateError(){
setLoading('errorSpinner',true);
const subst=parseFloat(document.getElementById('substRate').value)/100;
const del=parseFloat(document.getElementById('delRate').value)/100;
const loss=parseFloat(document.getElementById('lossRate').value)/100;
try{
const res=await fetch('/error',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({subst,del,loss})});
const data=await res.json();
if(data.error){showAlert('errorAlert','❌ '+data.error,'error');setLoading('errorSpinner',false);return}
if(data.success)showAlert('errorAlert','✅ Phục hồi hoàn toàn!','success');
else if(data.partial)showAlert('errorAlert','⚠️ Phục hồi một phần','warning');
else showAlert('errorAlert','❌ Không phục hồi được','error');
document.getElementById('errorMetrics').style.display='grid';
document.getElementById('errorMetrics').innerHTML=`<div class="metric"><div class="metric-value">${data.valid}</div><div class="metric-label">Hợp lệ</div></div><div class="metric"><div class="metric-value">${data.bad}</div><div class="metric-label">Hỏng/mất</div></div>`;
}catch(e){showAlert('errorAlert','❌ Lỗi','error')}
finally{setLoading('errorSpinner',false)}
}

async function runComparison(){
const loss=parseFloat(document.getElementById('compareLossRate').value)/100;
setLoading('compareSpinner',true);
try{
const res=await fetch('/compare',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({loss})});
const data=await res.json();
const ctx=document.getElementById('compareChart').getContext('2d');
if(chartInstance)chartInstance.destroy();
chartInstance=new Chart(ctx,{type:'bar',data:{labels:['Không sửa','Mã lặp','Fountain'],datasets:[{label:'Phục hồi (%)',data:[data.no,data.rep,data.fountain],backgroundColor:['#ff4d6d','#ffd166','#00ff9d']}]},options:{scales:{y:{max:100}}}})
}catch(e){}
finally{setLoading('compareSpinner',false)}
}

async function simulateReality(){
const enzyme=document.getElementById('enzymeSelect').value;
const temp=parseFloat(document.getElementById('tempSlider').value);
setLoading('realitySpinner',true);
try{
const res=await fetch('/reality',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({enzyme,temperature:temp})});
const data=await res.json();
if(data.error){showAlert('realityAlert','❌ '+data.error,'error');setLoading('realitySpinner',false);return}
if(data.success)showAlert('realityAlert','✅ Phục hồi hoàn toàn!','success');
else if(data.partial)showAlert('realityAlert','⚠️ Phục hồi một phần','warning');
else showAlert('realityAlert','❌ Không phục hồi được','error');
document.getElementById('realityMetrics').style.display='grid';
document.getElementById('realityMetrics').innerHTML=`<div class="metric"><div class="metric-value">${data.valid}</div><div class="metric-label">Hợp lệ</div></div><div class="metric"><div class="metric-value">${data.bad}</div><div class="metric-label">Hỏng/mất</div></div>`;
}catch(e){showAlert('realityAlert','❌ Lỗi','error')}
finally{setLoading('realitySpinner',false)}
}

async function runPCR(){
const cycles=parseInt(document.getElementById('pcrCycles').value);
setLoading('pcrSpinner',true);
try{
const res=await fetch('/pcr',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({cycles})});
const data=await res.json();
if(data.error){showAlert('pcrAlert','❌ '+data.error,'error');setLoading('pcrSpinner',false);return}
showAlert('pcrAlert','✅ '+data.total_copies+' bản sao, '+data.total_errors+' lỗi','success');
const ctx=document.getElementById('pcrChart').getContext('2d');
if(pcrChartInstance)pcrChartInstance.destroy();
pcrChartInstance=new Chart(ctx,{type:'line',data:{labels:data.labels,datasets:[{label:'Số bản sao',data:data.copies,borderColor:'#00ff9d',fill:false}]},options:{scales:{y:{type:'logarithmic'}}}})
}catch(e){showAlert('pcrAlert','❌ Lỗi','error')}
finally{setLoading('pcrSpinner',false)}
}

// 3D
function init3DDNA(){
if(!window._dna3dData) return;
const container=document.getElementById('dna3dContainer');
container.innerHTML='';
const w=container.clientWidth, h=300;
threeScene=new THREE.Scene();
threeCamera=new THREE.PerspectiveCamera(60,w/h,0.1,100);
threeCamera.position.set(3,2,8);
threeRenderer=new THREE.WebGLRenderer({antialias:true});
threeRenderer.setSize(w,h);
threeRenderer.setClearColor(0x1a2235);
container.appendChild(threeRenderer.domElement);
const light=new THREE.DirectionalLight(0xffffff,1);
light.position.set(5,5,5);
threeScene.add(light);
threeScene.add(new THREE.AmbientLight(0x404040));
const colors={A:0x00ff9d,T:0xff4d6d,C:0x00bfff,G:0xffd166};
window._dna3dData.forEach(p=>{
const geo=new THREE.SphereGeometry(0.15,16,16);
const mat=new THREE.MeshPhongMaterial({color:colors[p.nucleotide]||0xffffff});
const sphere=new THREE.Mesh(geo,mat);
sphere.position.set(p.x,p.y+2,p.z);
threeScene.add(sphere);
});
function animate(){requestAnimationFrame(animate);threeScene.rotation.y+=0.002;threeRenderer.render(threeScene,threeCamera)}
animate();
let isDragging=false, prevX, prevY;
container.addEventListener('touchstart',e=>{isDragging=true;prevX=e.touches[0].clientX;prevY=e.touches[0].clientY});
container.addEventListener('touchmove',e=>{if(!isDragging)return;const dx=e.touches[0].clientX-prevX,dy=e.touches[0].clientY-prevY;threeScene.rotation.y+=dx*0.01;threeScene.rotation.x+=dy*0.01;prevX=e.touches[0].clientX;prevY=e.touches[0].clientY});
container.addEventListener('touchend',()=>{isDragging=false});
}
document.addEventListener('DOMContentLoaded',()=>{
const observer=new MutationObserver(()=>{const tab3d=document.getElementById('tab-3d');if(tab3d&&tab3d.classList.contains('active')&&window._dna3dData)init3DDNA()});
document.querySelectorAll('.panel').forEach(p=>observer.observe(p,{attributes:true,attributeFilter:['class']}));
});

async function loadHistory(){
try{const res=await fetch('/history');const data=await res.json();if(data.history.length){let h='<table class="history-table"><tr><th>Thời gian</th><th>File</th><th>Kích thước</th><th>Droplet</th><th>GC</th></tr>';data.history.forEach(r=>h+=`<tr><td>${r.time}</td><td>${r.filename}</td><td>${r.size}</td><td>${r.droplets}</td><td>${r.gc}%</td></tr>`);h+='</table>';document.getElementById('historyContent').innerHTML=h}else document.getElementById('historyContent').innerHTML='<p style="color:var(--muted)">Chưa có lịch sử</p>'}catch(e){}
}

const devices=[{n:'DNA (lý tưởng)',y:1000},{n:'Băng từ',y:30},{n:'USB',y:10},{n:'SSD',y:5},{n:'HDD',y:5}];
document.getElementById('lifeChart').innerHTML=devices.map(d=>`<div style="display:flex;align-items:center;gap:.5rem;margin:.5rem 0"><span style="width:120px;font-size:.8rem">${d.n}</span><div style="flex:1;background:var(--surface2);height:20px;border-radius:4px;overflow:hidden"><div style="height:100%;width:${Math.max(d.y/10,3)}%;background:${d.y===1000?'#00ff9d':'#ff4d6d'};border-radius:4px;display:flex;align-items:center;justify-content:flex-end;padding-right:5px;font-size:.7rem;color:#0a0e1a">${d.y} năm</div></div></div>`).join('');
</script>
</body>
</html>"""

# ---------- ROUTES ----------
@app.route('/')
def index(): return HTML

@app.route('/encode_text', methods=['POST'])
def encode_text():
    try:
        data = request.json
        text = data.get('text','')
        opt = data.get('opt_mode','ai')
        huffman = data.get('huffman',True)
        pwd = data.get('password','')
        raw = text.encode()
        return _encode_common(raw, 'vanban.txt', opt, huffman, pwd)
    except Exception as e:
        return jsonify({'error': str(e)})

@app.route('/encode_file', methods=['POST'])
def encode_file():
    try:
        f = request.files.get('file')
        if not f: return jsonify({'error':'Không có file'})
        raw = f.read()
        opt = request.form.get('opt_mode','ai')
        huffman = request.form.get('huffman','true')=='true'
        pwd = request.form.get('password','')
        return _encode_common(raw, f.filename, opt, huffman, pwd)
    except Exception as e:
        return jsonify({'error': str(e)})

def _encode_common(raw, filename, opt_mode, use_huffman, password):
    stored['original_data'] = raw
    stored['password'] = password
    stored['use_huffman'] = use_huffman
    if password:
        raw = aes_encrypt(raw, password)
    if use_huffman:
        bits, codes, orig_len = huffman_compress(raw)
        padded = bits + '0'*((8-len(bits)%8)%8)
        compressed = bits_to_bytes(padded)
        stored['huffman_codes'] = codes
        stored['huffman_orig_len'] = orig_len
        data_to_store = compressed
    else:
        data_to_store = raw
    system = DNAStorageSystem(chunk_size=4, droplet_factor=3.0, scramble_trials=10)
    meta = system.store(data_to_store, filename)
    stored['meta'] = meta
    stored['data'] = data_to_store
    # 3D
    first_seq = meta['dna_sequences'][0] if meta['dna_sequences'] else ''
    dna3d = generate_dna_3d_coordinates(first_seq[:100])
    # Lịch sử
    history.append({'time':datetime.now().strftime('%H:%M:%S %d/%m/%Y'),'filename':filename,'size':len(raw),'droplets':meta['num_droplets'],'gc':round(meta['gc_content_avg']*100,1)})
    if len(history)>50: history.pop(0)
    preview = '\n'.join([f"Droplet {i}: {seq[:50]}... (dài {len(seq)})" for i,seq in enumerate(meta['dna_sequences'][:5])])
    return jsonify({
        'num_droplets': meta['num_droplets'],
        'total_nucleotides': meta['total_nucleotides'],
        'gc_content_avg': meta['gc_content_avg'],
        'max_homopolymer_max': meta['max_homopolymer_max'],
        'dna_preview': preview,
        'min_droplets': meta['num_blocks'],
        'dna_3d': dna3d
    })

@app.route('/decode', methods=['POST'])
def decode():
    if 'meta' not in stored: return jsonify({'error':'Chưa mã hoá'})
    pwd = request.json.get('password','')
    meta = stored['meta']
    system = DNAStorageSystem(chunk_size=meta['chunk_size'], droplet_factor=3.0)
    rec, _ = system.retrieve(meta['dna_sequences'], meta)
    if rec is None: return jsonify({'error':'Không đủ dữ liệu'})
    if stored.get('use_huffman') and stored.get('huffman_codes'):
        bits = bytes_to_bits(rec)
        rec = huffman_decompress(bits, stored['huffman_codes'], stored['huffman_orig_len'])
    if stored.get('password'):
        if pwd != stored['password']: return jsonify({'error':'Sai mật khẩu'})
        rec = aes_decrypt(rec, pwd)
    stored['recovered'] = rec
    try: text = rec.decode()
    except: text = None
    return jsonify({'text': text})

@app.route('/download')
def download():
    if 'recovered' not in stored: return 'Không có dữ liệu',404
    meta = stored.get('meta',{})
    filename = meta.get('filename','output.bin')
    return send_file(io.BytesIO(stored['recovered']), as_attachment=True, download_name=filename)

@app.route('/error', methods=['POST'])
def error_sim():
    if 'meta' not in stored: return jsonify({'error':'Chưa mã hoá'})
    body = request.json
    subst, dell, loss = body.get('subst',0), body.get('del',0), body.get('loss',0)
    meta = stored['meta']
    corrupted = [simulate_errors(s, subst, dell) for s in meta['dna_sequences']]
    system = DNAStorageSystem(chunk_size=meta['chunk_size'], droplet_factor=3.0)
    rec, stats = system.retrieve(corrupted, meta, error_loss=loss)
    result = {'valid':stats['valid'], 'bad':stats['corrupted']+stats['lost'], 'success':False, 'partial':False}
    if rec is not None:
        if rec == stored['data']: result['success'] = True
        else: result['partial'] = True
    return jsonify(result)

@app.route('/compare', methods=['POST'])
def compare():
    if 'data' not in stored: return jsonify({'error':'Chưa mã hoá'})
    loss = request.json.get('loss',0.1)
    data = stored['data']
    # No correction
    noc = 100 if random.random() >= loss else 0
    # Repetition
    rep = 100 if random.random() >= loss**3 else 0
    # Fountain
    system = DNAStorageSystem(chunk_size=4, droplet_factor=3.0)
    meta = system.store(data,'tmp')
    corrupted = [simulate_errors(s,0,0) for s in meta['dna_sequences']]
    rec, _ = system.retrieve(corrupted, meta, error_loss=loss)
    fountain = 100 if rec == data else 0
    return jsonify({'no':noc, 'rep':rep, 'fountain':fountain})

@app.route('/reality', methods=['POST'])
def reality():
    if 'meta' not in stored: return jsonify({'error':'Chưa mã hoá'})
    body = request.json
    enzyme = body.get('enzyme','Phusion')
    temp = body.get('temperature',55)
    rates = get_error_rates_from_enzyme(enzyme, temp)
    meta = stored['meta']
    corrupted = [simulate_errors(s, rates['subst'], rates['del'], rates['ins']) for s in meta['dna_sequences']]
    system = DNAStorageSystem(chunk_size=meta['chunk_size'], droplet_factor=3.0)
    rec, stats = system.retrieve(corrupted, meta, error_loss=0)
    result = {'valid':stats['valid'], 'bad':stats['corrupted']+stats['lost'], 'success':False, 'partial':False}
    if rec is not None:
        if rec == stored['data']: result['success'] = True
        else: result['partial'] = True
    return jsonify(result)

@app.route('/pcr', methods=['POST'])
def pcr():
    if 'meta' not in stored: return jsonify({'error':'Chưa mã hoá'})
    cycles = request.json.get('cycles',3)
    seq = stored['meta']['dna_sequences'][0]
    copies, errs = simulate_pcr(seq, cycles, 0.001)
    labels = [f'CK {i}' for i in range(cycles+1)]
    copy_counts = [1]
    for i in range(cycles): copy_counts.append(copy_counts[-1]*2)
    return jsonify({'total_copies':len(copies), 'total_errors':errs, 'labels':labels, 'copies':copy_counts})

@app.route('/history')
def hist(): return jsonify({'history':history})

# Import bổ sung cho decode
from dna_utils import bytes_to_bits

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8501, debug=False)