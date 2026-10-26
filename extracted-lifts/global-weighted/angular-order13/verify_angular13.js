"use strict";
// Exact coefficient checks; no dependencies. Run: node verify_angular13.js
function runAudit() {
 const add=(p,a,b,c)=>{let k=a+","+b;p.set(k,(p.get(k)||0n)+BigInt(c));if(!p.get(k))p.delete(k);};
 const plus=(p,q)=>{for(const [k,c]of q){let[a,b]=k.split(",").map(Number);add(p,a,b,c);}return p;};
 const square=z=>{let p=new Map();for(const[a,b,c]of z.t)for(const[d,e,f]of z.t)add(p,a+d+z.pu,b+e+z.pt,z.w*c*f);return p;};
 const eq=(p,q)=>p.size===q.size&&[...p].every(([k,v])=>q.get(k)===v);
 const poly=ts=>{let p=new Map();for(const[a,b,c]of ts)add(p,a,b,c);return p;};
 const ft=[[[0,0],[0,1],[0,2],[1,0],[1,1],[1,2]],[[0,0],[0,1],[0,2]],[[0,0],[0,1]]];
 const root=(w,pu,pt,t)=>({w,pu,pt,t});
 const b3={name:"triangle111",target:[[2,4,1],[1,0,1],[0,2,1],[1,2,-3]],roots:[root(1,0,0,[[1,2,1],[0,1,-1]]),root(1,1,0,[[0,1,1],[0,0,-1]]),root(2,1,1,[[0,1,1],[0,0,-1]])]};
 const bv={name:"triangle112_vertical",target:[[2,4,1],[0,4,1],[1,0,2],[1,2,-4]],roots:[root(1,0,0,[[1,2,1],[0,2,-1]]),root(2,1,0,[[0,2,1],[0,0,-1]])]};
 const bh={name:"triangle112_horizontal_repaired",target:[[2,4,1],[2,0,1],[0,2,2],[1,2,-4]],roots:[root(1,0,0,[[1,0,1],[0,1,-1]]),root(1,0,0,[[1,2,1],[0,1,-1]]),root(2,1,1,[[0,1,1],[0,0,-1]])]};
 function transform(B,s,r){
  let roots=B.roots.map(z=>{
   let t=z.t.map(([a,b,c])=>[a*s,b*r,c]),ma=Math.min(...t.map(x=>x[0])),mb=Math.min(...t.map(x=>x[1]));
   let pu=1-s+z.pu*s+2*ma,pt=2-2*r+z.pt*r+2*mb;
   let ua=Math.floor(pu/2),tb=Math.floor(pt/2);pu-=2*ua;pt-=2*tb;
   t=t.map(([a,b,c])=>[a-ma+ua,b-mb+tb,c]);return root(z.w,pu,pt,t);
  });
  return {name:B.name+"_"+s+"_"+r,roots,target:B.target.map(([a,b,c])=>[1-s+a*s,2-2*r+b*r,c])};
 }
 let cases=[];
 for(let s of[-1,1])for(let r of[-1,1])cases.push(transform(b3,s,r));
 for(let r of[-1,1])cases.push(transform(bv,1,r));
 for(let s of[-1,1])cases.push(transform(bh,s,1));
 for(let[a,b]of[[0,1],[1,0],[1,1],[1,-1]])cases.push({name:"pair_"+a+"_"+b,roots:[root(1,a===0?1:0,0,a===0?[[0,1+b,1],[0,1-b,-1]]:[[(1+a)/2,1+b,1],[(1-a)/2,1-b,-1]])],target:[[1+a,2+2*b,1],[1-a,2-2*b,1],[1,2,-2]]});
 const zero=n=>Array.from({length:n},()=>Array(n).fill(0));
 const gn=["g00","g01","g02","g10","h","g12","g20","g21","g22"];
 const an=[...Array.from({length:3},(_,i)=>Array.from({length:3},(_,j)=>"J"+i+j)).flat(),"C00","C01","C11","rho"];
 function pencil(x){
  let v=n=>x[n]||0,M=zero(6),B=zero(3),C=zero(2);
  ["g00","g01","g02","g20","g21","g22"].forEach((n,i)=>M[i][i]=2*v(n));
  for(let i=0;i<3;i++)for(let j=0;j<3;j++)M[i][3+j]=M[3+j][i]=2*v("J"+i+j);
  B[0][0]=2*v("g10")-4*v("J00");B[0][1]=B[1][0]=-2*v("J01")-2*v("J10")-v("C00");
  B[0][2]=B[2][0]=2*v("rho");B[1][1]=2*v("h")-4*(v("J02")+v("J11")+v("J20"))-4*v("rho")-4*v("C01");
  B[1][2]=B[2][1]=-2*v("J12")-2*v("J21")-v("C11");B[2][2]=2*v("g12")-4*v("J22");
  C[0][0]=2*v("C00");C[0][1]=C[1][0]=2*v("C01");C[1][1]=2*v("C11");
  return[M,B,C,[[2*v("g10")]],[[2*v("g12")]]];
 }
 function gramPoly(A,f,pu,pt){let p=new Map();for(let i=0;i<f.length;i++)for(let j=0;j<f.length;j++)add(p,f[i][0]+f[j][0]+pu,f[i][1]+f[j][1]+pt,A[i][j]);return p;}
 const evaluate=A=>plus(plus(gramPoly(A[0],ft[0],0,0),gramPoly(A[1],ft[1],1,0)),gramPoly(A[2],ft[2],1,1));
 let names=gn.concat(an),coo={};
 for(let n of names){
  let A=pencil({[n]:1}),p=new Map();
  if(gn.includes(n)){let k=gn.indexOf(n);add(p,Math.floor(k/3),2*(k%3),2);}
  if(!eq(evaluate(A),p))throw Error("affine identity: "+n);
  coo[n]=A.map(B=>B.flatMap((r,i)=>r.flatMap((v,j)=>i<=j&&v?[[i,j,v]]:[])));
 }
 const key=p=>[...p].filter(([k,v])=>v>0n).map(([k,v])=>{let[a,b]=k.split(",").map(Number);return[a-1,b/2-1,v.toString()].join(":");}).sort().join(";");
 let keys=new Set(),reports=[];
 for(let c of cases){
  let A=[zero(6),zero(3),zero(2)],sum=new Map();
  for(let z of c.roots){
   let bi=z.pu===0&&z.pt===0?0:z.pu===1&&z.pt===0?1:z.pu===1&&z.pt===1?2:-1;
   if(bi<0)throw Error("weight "+c.name);
   let f=ft[bi],v=Array(f.length).fill(0);
   for(let[a,b,k]of z.t){let i=f.findIndex(x=>x[0]===a&&x[1]===b);if(i<0)throw Error("span "+c.name+" "+a+","+b);v[i]+=k;}
   for(let i=0;i<v.length;i++)for(let j=0;j<v.length;j++)A[bi][i][j]+=z.w*v[i]*v[j];
   plus(sum,square(z));
  }
  let p=poly(c.target);if(!eq(sum,p))throw Error("circuit "+c.name);
  let x={};
  for(let[k,v]of p){let[a,b]=k.split(",").map(Number);if(b%2)throw Error("odd");x[gn[a*3+b/2]]=Number(v);}
  for(let i=0;i<3;i++)for(let j=0;j<3;j++)x["J"+i+j]=A[0][i][3+j];
  x.C00=A[2][0][0];x.C01=A[2][0][1];x.C11=A[2][1][1];x.rho=A[1][0][2];
  let P=pencil(x);
  if(!A.every((G,b)=>G.every((r,i)=>r.every((v,j)=>2*v===P[b][i][j]))))throw Error("diagonal restriction "+c.name);
  keys.add(key(p));reports.push({name:c.name,identity:true,fixed_span:true,fixed_diagonal_blocks:true,positive_outer_products:c.roots.length,certificate:x,roots:c.roots});
 }
 let E=[];for(let i=-1;i<=1;i++)for(let j=-1;j<=1;j++)if(i||j)E.push([i,j]);
 let det=(a,b)=>a[0]*b[1]-a[1]*b[0],gcd=(a,b)=>b?gcd(b,a%b):a,expected=new Set();
 for(let i=0;i<8;i++)for(let j=i+1;j<8;j++){
  if(E[i][0]+E[j][0]===0&&E[i][1]+E[j][1]===0){let p=new Map();for(let e of[E[i],E[j]])add(p,1+e[0],2+2*e[1],1);expected.add(key(p));}
  for(let k=j+1;k<8;k++){let w=[det(E[j],E[k]),det(E[k],E[i]),det(E[i],E[j])];if(w.every(v=>v<0))w=w.map(v=>-v);if(!w.every(v=>v>0))continue;let g=w.reduce(gcd);w=w.map(v=>v/g);let p=new Map();[i,j,k].forEach((n,l)=>add(p,1+E[n][0],2+2*E[n][1],w[l]));expected.add(key(p));}
 }
 if(expected.size!==12||keys.size!==12||[...expected].some(k=>!keys.has(k)))throw Error("enumeration");
 return{status:"passed",arithmetic:"exact integer coefficients (BigInt)",circuits_enumerated:12,circuit_checks:reports,affine_basis_identities:names.length,representation:{visible:gn,auxiliary:an,block_orders:[6,3,2,1,1],convention:"twice the affine pencil; blockwise upper-triangular [row,column,integer], zero based",coefficient_matrices:coo},proof_limit:"Coefficient identities and all circuit certificates checked. Global necessity and closure use the accompanying proof. Not a full Q222 affine lift."};
}
if(typeof process!=="undefined")console.log(JSON.stringify(runAudit(),null,2));
