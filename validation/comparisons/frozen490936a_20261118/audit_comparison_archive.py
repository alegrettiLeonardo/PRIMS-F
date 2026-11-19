#!/usr/bin/env python3
from pathlib import Path
import argparse,hashlib,json,zipfile
p=argparse.ArgumentParser();p.add_argument('--archive',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();pdfs=[];crcs={};hashes_verified=0
with zipfile.ZipFile(a.archive) as z:
 assert z.testzip() is None
 expected=json.loads(z.read('CONTENTS_SHA256.json'))['contents_sha256']
 for member in z.infolist():
  content=z.read(member);name=member.filename;crcs[name]=f'{member.CRC:08x}'
  assert not content.startswith((b'\x7fELF',b'MZ')),name
  assert Path(name).suffix.lower() not in {'.exe','.dll','.so','.o','.a','.mod'},name
  if name in expected:
   assert hashlib.sha256(content).hexdigest()==expected[name],name
   hashes_verified+=1
  if content.startswith(b'%PDF') or name.lower().endswith('.pdf'):
   assert content.startswith(b'%PDF') and b'Matplotlib' in content,name
   assert name.startswith('validation/comparisons/'),name
   assert Path(name).name in {'mean_cycle.pdf','osc_envelopes.pdf','osc_envelopes_last_cycles.pdf','osc_envelopes_spatial.pdf'} or '/stable_analysis/comparison_' in name,name
   pdfs.append(name)
 assert hashes_verified==len(expected)
result={'archive_sha256':hashlib.sha256(a.archive.read_bytes()).hexdigest(),'crc_test':'PASS',
 'member_count':len(crcs),'internal_file_sha256_verified':hashes_verified,'elf_pe_signature_scan':'PASS: none',
 'executable_extension_scan':'PASS: none','pdf_scope':'Only Matplotlib-generated comparison figures with allowed names; no bibliographic PDF',
 'generated_pdf_members':pdfs,'member_crc32':crcs}
a.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['member_crc32','generated_pdf_members']}))
