/* Exercise the generated report, including every linked scientific artifact. */
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import {JSDOM} from 'jsdom';

const root=path.resolve(process.argv[2]||'docs/vaccine-results');
const payload=JSON.parse(fs.readFileSync(path.join(root,'data.json'),'utf8'));
const dom=new JSDOM(fs.readFileSync(path.join(root,'index.html'),'utf8'),{url:'http://localhost/',runScripts:'outside-only'});
const w=dom.window;
w.fetch=async()=>({ok:true,json:async()=>payload});
w.HTMLElement.prototype.scrollIntoView=()=>{};
let copied='';Object.defineProperty(w.navigator,'clipboard',{value:{writeText:async s=>{copied=s;}}});
w.eval(fs.readFileSync(path.join(root,'app.js'),'utf8'));
await new Promise(resolve=>setImmediate(resolve));
const $=s=>w.document.querySelector(s), $$=s=>[...w.document.querySelectorAll(s)];
const input=(selector,value)=>{const e=$(selector);e.value=value;e.dispatchEvent(new w.Event('input',{bubbles:true}));};
assert.equal($('#load-error').hidden,true,$('#load-error').textContent);
for(const [key,d] of Object.entries(payload.designs)){
 $(`[data-mode="${key}"]`).click();
 input('#protein-search','');input('#hla-search','');$('[data-locus="all"]').click();
 assert.equal($$('#protein-rows tr').length,d.proteins.length);
 assert.equal($$('#hla-rows tr').length,d.panel_size);
 assert.match($('#summary-metrics').textContent,new RegExp(String(d.pmhc)));
 input('#coverage-prefix',0);assert.match($('#coverage-added').textContent,/No proteins/);
 input('#coverage-prefix',d.proteins.length);assert.match($('#coverage-added').textContent,/Latest addition/);
 const final=d.coverage.find(r=>r.axis==='protein'&&r.step===d.proteins.length);
 assert.equal(final.ms_peptides,d.ms_peptides);assert.equal(final.pmhc,d.pmhc);
 assert.equal(final.alleles,d.supported_alleles);
 assert.equal($$('#coverage-cancers tr').length,d.cancers.length);
 input('#protein-search','MAGE');assert.equal($$('#protein-rows tr').length,1);
 assert.match($('#protein-rows').textContent,/MAGEA4/);input('#protein-search','');
 input('#tissue-search','Heart');assert.ok($('#tissue-rows').textContent.includes('Heart')||$('#tissue-rows').textContent.includes('No matching'));
 input('#tissue-search','');
 $('[data-sequence="full"]').click();assert.equal($('#sequence-text').textContent.replace(/\s/g,''),d.sequences.full);
 $('#copy-sequence').click();await new Promise(resolve=>setImmediate(resolve));assert.equal(copied,d.sequences.full);
 for(const a of $$('#download-links a'))assert.ok(fs.existsSync(path.join(root,decodeURIComponent(a.getAttribute('href')))));
 for(const img of $$('img:not([hidden])')){
  assert.ok(fs.existsSync(path.join(root,img.getAttribute('src'))),img.getAttribute('src'));
  assert.ok(img.width>0&&img.height>0,'Scientific images must reserve their intrinsic aspect ratio');
 }
 assert.equal(d.sequences.protein.length,d.protein_aa);assert.equal(d.sequences.full.length,d.total_nt);
 assert.equal(d.hla.reduce((a,r)=>a+r.retained_peptides,0),d.pmhc);
 assert.equal(d.proteins.reduce((a,p)=>a+p.assembled_pieces,0),d.native_pieces);
 assert.ok(d.proteins.every(p=>p.name.split('/').every(g=>!g.startsWith('MAGE')||g==='MAGEA4')));
 assert.ok(d.protein_aa<=d.max_aa&&d.total_nt<=d.max_nt);
}
assert.doesNotMatch(w.document.body.textContent,/\b(PR|issues?|development|merged|deployed|PyPI)\b/i);
assert.match($('#ct83-comparison').textContent,/CT83|not screened|selected|budget/);
console.log(`Website verified: ${Object.keys(payload.designs).length} designs, coverage controls, source maps, counts, sequences and every downloadable artifact.`);
