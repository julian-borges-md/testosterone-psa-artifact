const fs=require('fs');
const {Document,Packer,Paragraph,TextRun,HeadingLevel,Table,TableRow,TableCell,WidthType,ShadingType,AlignmentType,BorderStyle,Footer,PageNumber,LevelFormat,ImageRun}=require('docx');
const sizeOf=(p)=>{const b=fs.readFileSync(p);let i=8;while(i<b.length){const len=b.readUInt32BE(i);const type=b.toString('ascii',i+4,i+8);if(type==='IHDR')return{w:b.readUInt32BE(i+8),h:b.readUInt32BE(i+12)};i+=12+len;}return{w:1200,h:700};};
const SRC=process.argv[2]||'manuscript.md';const OUT=process.argv[3]||'/mnt/user-data/outputs/Borges_Testosterone_Raises_PSA_Not_Cancer_Risk_v1.2.docx';const lines=fs.readFileSync(SRC,'utf8').split('\n');
const FONT='Times New Roman', SZ=24, TW=9360;
function runs(text,opt={}){ // **bold** segments
  const out=[]; const parts=text.split(/(\*\*[^*]+\*\*)/);
  for(const p of parts){ if(!p) continue;
    if(p.startsWith('**')) out.push(new TextRun({text:p.slice(2,-2),bold:true,font:FONT,size:opt.size||SZ}));
    else out.push(new TextRun({text:p,font:FONT,size:opt.size||SZ,italics:opt.italics})); }
  return out;}
const B={style:BorderStyle.SINGLE,size:4,color:'808080'};
function table(caption,rows){
  const cells=rows.map(r=>r.replace(/^\|\s*/,'').replace(/\s*\|$/,'').split(/\s*\|\s*/));
  const n=cells[0].length;
  const len=Array(n).fill(0); cells.forEach(r=>r.forEach((c,i)=>len[i]=Math.max(len[i],Math.min(Math.max(c.length,13),40))));
  const tot=len.reduce((a,b)=>a+b,0); let w=len.map(l=>Math.max(1250,Math.round(TW*l/tot)));
  if(caption.startsWith('Table 1.')) w=[1250,2550,2400,1550,1610];
  const s=w.reduce((a,b)=>a+b,0); w=w.map(x=>Math.floor(x*TW/s)); w[n-1]+=TW-w.reduce((a,b)=>a+b,0);
  const trs=cells.map((r,ri)=>new TableRow({tableHeader:ri===0,children:r.map((c,ci)=>new TableCell({
    width:{size:w[ci],type:WidthType.DXA},borders:{top:B,bottom:B,left:B,right:B},
    shading:ri===0?{fill:'E7E6E6',type:ShadingType.CLEAR,color:'auto'}:undefined,
    margins:{top:60,bottom:60,left:90,right:90},
    children:[new Paragraph({children:[new TextRun({text:c,bold:ri===0,font:FONT,size:18})]})]}))}));
  return [new Paragraph({spacing:{before:200,after:80},keepNext:true,children:[new TextRun({text:caption,bold:true,font:FONT,size:20})]}),
          new Table({width:{size:TW,type:WidthType.DXA},columnWidths:w,rows:trs}),
          new Paragraph({spacing:{after:120},children:[]})];
}
const kids=[]; let i=0;
while(i<lines.length){ const l=lines[i];
  if(!l.trim()){i++;continue;}
  if(l.startsWith('# ')) kids.push(new Paragraph({heading:HeadingLevel.TITLE,alignment:AlignmentType.CENTER,spacing:{after:240},children:[new TextRun({text:l.slice(2),bold:true,font:FONT,size:32})]}));
  else if(l.startsWith('### ')) kids.push(new Paragraph({heading:HeadingLevel.HEADING_2,spacing:{before:200,after:80},keepNext:true,children:[new TextRun({text:l.slice(4),bold:true,font:FONT,size:24,color:'000000'})]}));
  else if(l.startsWith('## ')) kids.push(new Paragraph({heading:HeadingLevel.HEADING_1,spacing:{before:300,after:120},keepNext:true,children:[new TextRun({text:l.slice(3),bold:true,font:FONT,size:28,color:'000000'})]}));
  else if(l.startsWith('TABLE:')){ const cap=l.slice(6); const rows=[]; i++; while(lines[i]!=='END'){rows.push(lines[i]);i++;} kids.push(...table(cap,rows)); }
  else if(l.startsWith('EQ:')) kids.push(new Paragraph({alignment:AlignmentType.CENTER,spacing:{before:120,after:120},children:[new TextRun({text:l.slice(3),italics:true,font:'Cambria Math',size:SZ})]}));
  else if(l.startsWith('FIG:')){ const [f,cap]=l.slice(4).split('|'); const {w,h}=sizeOf(f); const W=600, H=Math.round(600*h/w);
    kids.push(new Paragraph({alignment:AlignmentType.CENTER,spacing:{before:240,after:80},children:[new ImageRun({data:fs.readFileSync(f),transformation:{width:W,height:H}})]}));
    kids.push(new Paragraph({spacing:{after:240},children:[new TextRun({text:cap,font:FONT,size:18})]})); }
  else if(l.startsWith('REF:')) kids.push(new Paragraph({numbering:{reference:'refs',level:0},spacing:{after:80},children:runs(l.slice(4),{size:22})}));
  else kids.push(new Paragraph({spacing:{after:160,line:300},alignment:AlignmentType.JUSTIFIED,children:runs(l)}));
  i++; }
const doc=new Document({creator:'Julian Yin Vieira Borges',title:'Testosterone raises PSA, not prostate cancer risk',
  styles:{default:{document:{run:{font:FONT,size:SZ}}}},
  numbering:{config:[{reference:'refs',levels:[{level:0,format:LevelFormat.DECIMAL,text:'%1.',alignment:AlignmentType.LEFT,style:{paragraph:{indent:{left:540,hanging:540}}}}]}]},
  sections:[{properties:{page:{size:{width:12240,height:15840},margin:{top:1440,bottom:1440,left:1440,right:1440}}},
    footers:{default:new Footer({children:[new Paragraph({alignment:AlignmentType.CENTER,children:[new TextRun({children:[PageNumber.CURRENT],font:FONT,size:20})]})]})},
    children:kids}]});
Packer.toBuffer(doc).then(b=>{fs.writeFileSync(OUT,b);console.log('ok',b.length);});
