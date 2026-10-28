#!/usr/bin/env python3
"""Independent exact replay of the small residual group tail only.

This does not import either group producer or assert that the 60,485-move
prefix has been independently replayed here. The parent audit owns that step.
"""
from pathlib import Path
from collections import Counter
import copy,gzip,hashlib,json

ROOT=Path(__file__).resolve().parent
GROUP=ROOT.parent/'group_audit'
def require(ok,msg):
    if not ok:raise ValueError(msg)
def load(path):return json.loads(gzip.decompress(path.read_bytes()))
def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def inverse(word):return tuple(-v for v in word[::-1])
def reduced(word):
    stack=[]
    for letter in word:
        if stack and stack[-1]+letter==0:stack.pop()
        else:stack.append(letter)
    first=0;last=len(stack)
    while last-first>1 and stack[first]+stack[last-1]==0:first+=1;last-=1
    return tuple(stack[first:last])
def rotate(word,n):return word[n:]+word[:n]

def replay(prefix,tail):
    require(tail['initial_generators']==prefix['generators'],'tail starting generators')
    require(tail['initial_relators_sha256']==digest(prefix['relators']),'tail starting relators hash')
    rels=set(map(tuple,prefix['relators']));generators=set(prefix['generators']);counts=Counter()
    for event in tail['events']:
        kind=event['kind'];counts[kind]+=1
        if kind=='relator_shortening':
            rule=tuple(event['rule_relator']);target=tuple(event['target'])
            require(rule in rels and target in rels and rule!=target,'shortening rule must be separately retained')
            sign=event['rule_sign'];require(sign in (-1,1),'rule sign')
            j=event['rule_start'];p=event['match_start'];n=event['match_length']
            require(0<=j<len(rule) and 0<=p<len(target) and 0<n<=min(len(rule),len(target)),'shortening indices')
            relation=rotate(rule if sign==1 else inverse(rule),j);word=rotate(target,p)
            require(relation[:n]==word[:n],'shortening does not match the retained rule')
            expected=reduced(inverse(relation[n:])+word[n:])
            require(expected==tuple(event['result']) and len(expected)<len(target),'shortening result')
            rels.remove(target)
            if expected:rels.add(expected)
        elif kind=='eliminate':
            g=event['generator'];rule=tuple(event['relator']);replacement=tuple(event['replacement'])
            require(g in generators and rule in rels and sum(abs(v)==g for v in rule)==1,'defining generator relation')
            i=next(i for i,v in enumerate(rule) if abs(v)==g);rest=rule[i+1:]+rule[:i]
            require(replacement==(inverse(rest) if rule[i]>0 else rest),'generator replacement')
            rels.remove(rule);opposite=inverse(replacement)
            rels={reduced(v for letter in word for v in (replacement if letter==g else opposite if letter==-g else (letter,))) for word in rels}
            rels.discard(());generators.remove(g)
        else:raise ValueError('unknown tail move')
        require(all(abs(v) in generators for w in rels for v in w),'eliminated generator remains')
    require(sorted(generators)==tail['generators'] and sorted(rels)==list(map(tuple,tail['relators'])),'tail endpoint')
    require(not generators and not rels,'tail did not reach the empty presentation')
    return dict(counts)

if __name__=='__main__':
    prefix_path=GROUP/'tietze_certificate.json.gz';tail_path=GROUP/'residual_tietze_certificate.json.gz'
    prefix=load(prefix_path);tail=load(tail_path)
    require(tail['prefix_file_sha256']==hashlib.sha256(prefix_path.read_bytes()).hexdigest(),'tail prefix bytes')
    counts=replay(prefix,tail)
    a,b=prefix['generators'];words={47:(-b,a,a,b,-a),46:(-b,a,a,a,b,-a,-a,-a),262:(b,b,-a,-b,-b,-b,a,b),26:(-b,a,-b,a,-b,-a,b,b,-a)}
    for i,w in words.items():require(tuple(prefix['relators'][i])==w,'short algebra source relation')
    wrong=copy.deepcopy(tail);wrong['events'][5]['result'][0]*=-1
    try:replay(prefix,wrong)
    except ValueError as e:negative={'name':'corrupt_shortened_relator','rejected':True,'reason':str(e)}
    else:raise ValueError('corrupt tail escaped')
    receipt={'schema':'independent-residual-group-tail-review-v1','result':'PASS','scope':'Residual two-generator tail only; source presentation and 60,485-move prefix replay belong to the parent group audit.','prefix_file_sha256':hashlib.sha256(prefix_path.read_bytes()).hexdigest(),'tail_file_sha256':hashlib.sha256(tail_path.read_bytes()).hexdigest(),'initial_generators':prefix['generators'],'initial_relators':len(prefix['relators']),'moves':counts,'endpoint':'empty presentation','short_algebra_original_relator_indices':list(words),'negative_control':negative}
    (ROOT/'group_tail_review.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
