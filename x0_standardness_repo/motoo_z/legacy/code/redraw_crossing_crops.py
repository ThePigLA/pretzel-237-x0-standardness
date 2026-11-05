#!/usr/bin/env python3
"""Source-PDF crops for manual review of a guided crossing record."""
import argparse, json, math
from pathlib import Path
import fitz
from PIL import Image, ImageDraw, ImageFont

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('record',type=Path);p.add_argument('--pdf',required=True,type=Path)
p.add_argument('--out',required=True,type=Path);p.add_argument('--ids',nargs='*',type=int)
p.add_argument('--all',action='store_true');p.add_argument('--radius',type=float,default=23)
p.add_argument('--scale',type=float,default=5)
a=p.parse_args();record=json.loads(a.record.read_text());page=fitz.open(a.pdf)[0]
ids=set(a.ids if a.ids is not None else (range(len(record['crossings'])) if a.all else record['unresolved_crossings']))
rows=[c for c in record['crossings'] if c['id'] in ids];a.out.mkdir(parents=True,exist_ok=True)
font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',14)
tiles=[];size=int(2*a.radius*a.scale);header=50
for c in rows:
 x,y=c['position_pdf'];clip=fitz.Rect(x-a.radius,y-a.radius,x+a.radius,y+a.radius)
 pix=page.get_pixmap(matrix=fitz.Matrix(a.scale,a.scale),clip=clip,alpha=False)
 im=Image.frombytes('RGB',(pix.width,pix.height),pix.samples)
 tile=Image.new('RGB',(size,header+size),'white');tile.paste(im,(0,header));d=ImageDraw.Draw(tile)
 names='/'.join(v['colour'] for v in c['visits'])
 d.text((4,2),f"{c['id']}: {names}",fill='black',font=font)
 d.text((4,20),f"({x:.1f}, {y:.1f}) over={c['over_visit']}",fill='black',font=font)
 # Crossing is centred. Corner guides leave the crossing ink untouched.
 cx,cy=size/2,header+size/2
 for dx,dy in [(-1,-1),(-1,1),(1,-1),(1,1)]:
  d.line((cx+dx*10,cy+dy*10,cx+dx*18,cy+dy*18),fill=(120,120,120),width=1)
 tile.save(a.out/f"crossing_{c['id']:03d}.png");tiles.append(tile)
for first in range(0,len(tiles),20):
 batch=tiles[first:first+20];cols=min(4,len(batch));nrows=math.ceil(len(batch)/cols)
 sheet=Image.new('RGB',(cols*(size+8),nrows*(size+header+8)),(220,220,220))
 for n,tile in enumerate(batch):sheet.paste(tile,((n%cols)*(size+8),(n//cols)*(size+header+8)))
 sheet.save(a.out/f'sheet_{first//20:02d}.png')
print(json.dumps({'crossing_crops':len(tiles),'output':str(a.out)}))
