#!/usr/bin/env python3
"""Portable macOS/Linux wall-time bound with live log and an attempt receipt.

This wrapper does not certify native computations. It preserves the command's
completed output files and records whether the process returned or timed out.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds',type=float,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('command',nargs=argparse.REMAINDER)
    args=parser.parse_args()
    command=args.command[1:] if args.command[:1]==['--'] else args.command
    if not command or args.seconds<=0:
        parser.error('provide a positive bound and a command after --')
    args.out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic()
    receipt={'command':command,'cwd':os.getcwd(),'wall_seconds_limit':args.seconds,
             'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'certificate_asserted':False}
    process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                             text=True,bufsize=1,start_new_session=True)
    def copy_output():
        with (args.out/'process.log').open('w') as stream:
            for line in process.stdout:
                stream.write(line);stream.flush()
                sys.stdout.write(line);sys.stdout.flush()
    reader=threading.Thread(target=copy_output)
    reader.start()
    timed_out=False
    interrupted=False
    try:
        process.wait(timeout=args.seconds)
    except (subprocess.TimeoutExpired,KeyboardInterrupt) as exc:
        timed_out=isinstance(exc,subprocess.TimeoutExpired)
        interrupted=not timed_out
        os.killpg(process.pid,signal.SIGTERM)
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGKILL)
            process.wait()
    reader.join()
    receipt.update({'elapsed_seconds':time.monotonic()-start,'returncode':process.returncode,
                    'hard_timeout':timed_out,'interrupted':interrupted,
                    'status':'TIMED_OUT_CHECK_SAVED_CHECKPOINTS' if timed_out else
                             'INTERRUPTED_CHECK_SAVED_CHECKPOINTS' if interrupted else
                             'COMMAND_RETURNED'})
    (args.out/'ATTEMPT_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return 124 if timed_out else 130 if interrupted else process.returncode


if __name__=='__main__':
    raise SystemExit(main())
