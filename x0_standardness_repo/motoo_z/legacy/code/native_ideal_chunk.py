#!/usr/bin/env python3
"""Bounded native ideal-body simplification with a saved returned checkpoint.

Only Regina's atomic progress tracker is read while its worker mutates the
triangulation. Individual moves are not logged; this is a search/checkpoint
lane, not an explicit-move certificate.
"""
import argparse,hashlib,json,threading,time
from pathlib import Path
import regina

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('input',type=Path);p.add_argument('--out',type=Path,required=True)
p.add_argument('--seconds',type=float,default=240);p.add_argument('--seed',type=int,default=42)
p.add_argument('--mode',choices=['full','greedy','updown'],default='full')
p.add_argument('--max24',type=int,default=4);p.add_argument('--max33',type=int,default=20)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
tri=regina.open(str(a.input));
if not isinstance(tri,regina.Triangulation4):raise ValueError('Need a 4D triangulation')
tri=regina.Triangulation4(tri);regina.RandomEngine.reseedWithDefault()
if hasattr(regina.RandomEngine,'reseed'):regina.RandomEngine.reseed(a.seed)
before=tri.size();record={'source_file':str(a.input.resolve()),'source_sha256':sha(a.input),
 'source_pentachora':before,'mode':a.mode,'seconds':a.seconds,'seed':a.seed,
 'working_interpretation':'Figure2 red over green dashed back arc; user adopted this assumption on 2026-11-02; author confirmation pending.',
 'individual_move_log_available':False,'standardness_proved':False,'closed_completion_asserted':False}
save(a.out/'START.json',record)
tracker=regina.ProgressTrackerObjective(before);answer={};start=time.monotonic()
def work():
 try:
  if a.mode=='full':answer['changed']=bool(tri.simplify(tracker))
  elif a.mode=='greedy':answer['changed']=bool(tri.simplifyToLocalMinimum())
  else:answer['changed']=bool(tri.simplifyUpDown(a.max24,a.max33,False))
 except BaseException as exc:answer['error']=repr(exc)
w=threading.Thread(target=work);w.start();last=start
with (a.out/'progress.jsonl').open('w') as log:
 while w.is_alive():
  w.join(1);elapsed=time.monotonic()-start
  if a.mode=='full' and elapsed>=a.seconds and not tracker.isCancelled():tracker.cancel()
  if time.monotonic()-last>=10:
   row={'elapsed_seconds':elapsed,'best_pentachora':tracker.objective() if a.mode=='full' else None,
        'stage':tracker.description() if a.mode=='full' else a.mode,'cancel_requested':tracker.isCancelled()}
   line=json.dumps(row);print(line,flush=True);log.write(line+'\n');log.flush();last=time.monotonic()
if 'error' in answer:raise RuntimeError(answer['error'])
packet=regina.PacketOfTriangulation4(tri);packet.setLabel('Figure2 working ideal native checkpoint')
if not packet.save(str(a.out/'checkpoint.rga')):raise RuntimeError('Save failed')
(a.out/'checkpoint.isosig').write_text(tri.isoSig()+'\n')
record.update(answer);record.update({'elapsed_seconds':time.monotonic()-start,'pentachora':tri.size(),
 'checkpoint_sha256':sha(a.out/'checkpoint.rga'),'isosig_sha256':sha(a.out/'checkpoint.isosig'),
 'cancel_requested':tracker.isCancelled(),'status':'NATIVE_CHECKPOINT_SAVED'})
save(a.out/'RESULT.json',record);print(json.dumps(record),flush=True)
