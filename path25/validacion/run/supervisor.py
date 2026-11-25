"""Exactly one hidden child; 60s/256MiB/32MiB ceilings; no retry."""
from pathlib import Path
import ctypes,json,subprocess,sys,time,hashlib
from ctypes import wintypes
D=Path(__file__).resolve().parent
design=json.loads((D/'inputs'/'DESIGN.json').read_text(encoding='utf-8'))
budget=design['resource_ceilings_for_first_run']
run=D/'attempt_01';run.mkdir(exist_ok=False)
class Memory(ctypes.Structure):
 _fields_=[('cb',wintypes.DWORD),('faults',wintypes.DWORD),('peak_rss',ctypes.c_size_t),
  ('rss',ctypes.c_size_t),('peak_paged',ctypes.c_size_t),('paged',ctypes.c_size_t),
  ('peak_nonpaged',ctypes.c_size_t),('nonpaged',ctypes.c_size_t),
  ('pagefile',ctypes.c_size_t),('peak_pagefile',ctypes.c_size_t)]
gm=ctypes.windll.psapi.GetProcessMemoryInfo
gm.argtypes=[wintypes.HANDLE,ctypes.POINTER(Memory),wintypes.DWORD];gm.restype=wintypes.BOOL
current=ctypes.windll.kernel32.GetCurrentProcess;current.restype=wintypes.HANDLE
def memory(handle):
 m=Memory();m.cb=ctypes.sizeof(m)
 if not gm(handle,ctypes.byref(m),m.cb):raise ctypes.WinError()
 return m.rss,m.peak_rss
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,x):(run/name).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
worker_hash=sha(D/'worker.py');supervisor_hash=sha(__file__)
save('GATE.json',{'worker_sha256':worker_hash,'supervisor_sha256':supervisor_hash,
 'design_sha256':sha(D/'inputs'/'DESIGN.json'),'attempt':1,'retries':0,'budgets':budget})
start=time.perf_counter();reason=None;peak_child=0;peak_parent=0;peak_sampled=0;last_disk=0;max_disk=0
with (run/'stdout.txt').open('wb') as stdout,(run/'stderr.txt').open('wb') as stderr:
 p=subprocess.Popen([sys.executable,'-B','-X','utf8',str(D/'worker.py'),str(run)],
    stdout=stdout,stderr=stderr,creationflags=subprocess.CREATE_NO_WINDOW)
 try:
  while p.poll() is None:
   now=time.perf_counter()
   try:
    cr,cp=memory(wintypes.HANDLE(int(p._handle)));pr,pp=memory(current())
    peak_child=max(peak_child,cp);peak_parent=max(peak_parent,pp);peak_sampled=max(peak_sampled,cr+pr)
   except OSError:
    if p.poll() is None:raise
   if now-start>budget['wall_seconds']:reason='WALL_LIMIT'
   if peak_child+peak_parent>budget['process_memory_MiB']*1024**2:reason='MEMORY_LIMIT'
   if now-last_disk>=.25:
    max_disk=max(max_disk,sum(x.stat().st_size for x in run.iterdir() if x.is_file()));last_disk=now
    if max_disk>budget['output_MiB']*1024**2:reason='OUTPUT_LIMIT'
   if reason:p.kill();break
   time.sleep(.05)
  p.wait(timeout=5)
  try:peak_child=max(peak_child,memory(wintypes.HANDLE(int(p._handle)))[1])
  except OSError:pass
 except BaseException:
  p.kill();p.wait();raise
elapsed=time.perf_counter()-start
disk=sum(x.stat().st_size for x in run.iterdir() if x.is_file())
if elapsed>budget['wall_seconds'] and reason is None:reason='FINAL_WALL_LIMIT'
if peak_child+peak_parent>budget['process_memory_MiB']*1024**2:reason='FINAL_MEMORY_LIMIT'
if disk>budget['output_MiB']*1024**2:reason='FINAL_OUTPUT_LIMIT'
status='PASS' if reason is None and p.returncode==0 and (run/'RESULTS.json').exists() else 'STOPPED_OR_FAILED'
report={'status':status,'returncode':p.returncode,'stop_reason':reason,'elapsed_seconds':elapsed,
 'peak_worker_rss_bytes':peak_child,'peak_supervisor_rss_bytes':peak_parent,
 'combined_peak_sum_upper_bound_bytes':peak_child+peak_parent,'peak_combined_sampled_bytes':peak_sampled,
 'output_bytes_before_supervision_receipt':disk,'budgets':budget,'attempts':1,'retries':0,
 'worker_sha256':worker_hash,'supervisor_sha256':supervisor_hash,
 'source_hashes_unchanged':sha(D/'worker.py')==worker_hash and sha(__file__)==supervisor_hash,
 'scope':'Execution resources only, not development time or performance comparison with other methods.'}
save('SUPERVISION.json',report);print(json.dumps(report),flush=True)
if status!='PASS':sys.exit(2)
