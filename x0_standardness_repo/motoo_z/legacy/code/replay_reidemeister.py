#!/usr/bin/env python3
"""Replay a marked Reidemeister certificate from its bound source diagram."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import regina
from normalize_diagram import metrics, perform, state


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('certificate',type=Path)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();raw=a.source.read_bytes();cert_raw=a.certificate.read_bytes()
    source=json.loads(raw);cert=json.loads(cert_raw)
    sha=hashlib.sha256(raw).hexdigest()
    if sha!=cert['source_json_sha256']:raise ValueError('Source hash mismatch')
    if source['decorations']!=cert['decorations']:raise ValueError('Handle decorations changed')
    order=source.get('pd_component_to_source',source.get('source_component_identity',list(range(len(source['decorations'])))))
    if order!=cert['source_component_identity']:raise ValueError('Marked component ordering changed')
    start=regina.Link.fromPD([[x+1 for x in row] for row in source['pd_code']])
    if state(start)!=cert['initial_state']:raise ValueError('Source diagram differs from certificate start')
    endpoint=start
    for i,record in enumerate(cert['moves']):
        if perform(endpoint,record['op'])!=record:raise ValueError(f'Move replay differs at {i}')
    if state(endpoint)!=cert['endpoint_state']:raise ValueError('Endpoint state differs')
    if 'regina_pd_code' in cert:
        if endpoint.pdData()!=cert['regina_pd_code']:raise ValueError('Exported endpoint PD differs')
    if metrics(endpoint,cert['decorations'])[1]!=cert['endpoint_metrics']:raise ValueError('Endpoint metrics differ')
    result={'status':'PASS_MARKED_REIDEMEISTER_REPLAY','source_sha256':sha,
            'certificate_sha256':hashlib.sha256(cert_raw).hexdigest(),
            'move_count':len(cert['moves']),'regina_distribution_version':importlib.metadata.version('regina'),
            'component_and_integer_framing_tracking':True,'endpoint_metrics':cert['endpoint_metrics'],
            'counterclockwise_planar_position_not_certified':True}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
