"use strict";
// Exact rational seed and selected-contact-minor checks. No dependencies.
function checkContactSeeds(){
 const gcd=(a,b)=>b?gcd(b,a%b):(a<0n?-a:a);
 const Q=(n,d=1n)=>{n=BigInt(n);d=BigInt(d);if(!d)throw Error("zero denominator");if(d<0n){n=-n;d=-d;}let g=gcd(n,d);return[n/g,d/g];};
 const add=(a,b)=>Q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),neg=a=>[-a[0],a[1]],sub=(a,b)=>add(a,neg(b));
 const mul=(a,b)=>Q(a[0]*b[0],a[1]*b[1]),div=(a,b)=>Q(a[0]*b[1],a[1]*b[0]),str=a=>a[1]===1n?""+a[0]:a[0]+"/"+a[1];
 const zero=()=>Q(0),one=()=>Q(1);
 const dual=(a,i=-1)=>Array.from({length:6},(_,j)=>j===0?a:j-1===i?one():zero());
 const da=(a,b)=>a.map((x,i)=>add(x,b[i])),dn=a=>a.map(neg),ds=(a,b)=>da(a,dn(b));
 const dm=(a,b)=>a.map((x,i)=>i?add(mul(a[i],b[0]),mul(a[0],b[i])):mul(a[0],b[0]));
 const dc=(n,d=1)=>dual(Q(n,d)),sq=a=>dm(a,a),scale=(n,a)=>dm(dc(n),a);
 function packets(X0,Y0){
  let W=dual(one(),0),X=dual(X0,1),Y=dual(Y0,2),u=dual(one(),3),v=dual(one(),4),Z=dc(1);
  let Du=da(dc(1),sq(u)),Dv=da(dc(1),sq(v));
  let Rm=ds(da(ds(dm(v,W),X),dm(dm(u,v),Y)),dm(u,Z));
  let Rd=ds(ds(da(dm(u,W),dm(dm(u,v),X)),Y),dm(v,Z));
  let al=dm(Du,sq(ds(dm(v,Y),Z))),be=dm(Dv,sq(da(X,dm(u,Z))));
  let core=scale(2,ds(da(da(sq(Rd),al),be),dm(dm(Du,Dv),sq(Z))));
  return[core,scale(4,dm(Du,sq(da(W,dm(v,X))))),scale(4,dm(dm(Du,Dv),sq(X))),scale(4,sq(Rm)),scale(4,dm(Dv,sq(ds(dm(u,W),Y)))),scale(4,dm(dm(Du,Dv),sq(Y))),scale(4,sq(Rd)),scale(4,al),scale(4,be)];
 }
 function det(M){
  let A=M.map(r=>r.slice()),d=one();
  for(let i=0;i<A.length;i++){
   let j=i;while(j<A.length&&!A[j][i][0])j++;if(j===A.length)return zero();
   if(j!==i){[A[j],A[i]]=[A[i],A[j]];d=neg(d);}
   let p=A[i][i];d=mul(d,p);
   for(let r=i+1;r<A.length;r++){let t=div(A[r][i],p);for(let c=i;c<A.length;c++)A[r][c]=sub(A[r][c],mul(t,A[i][c]));}
  }return d;
 }
 const seeds=[
  {name:"C1",X:Q(-1,2),Y:Q(1,4),c:[Q(27,2),Q(0),Q(0),Q(4),Q(2),Q(9),Q(5,4),Q(13,4),Q(9,4)],cols:[0,3,5,6,7,8],expected:"69120",expectedS:"103/37",expectedExpose:"107/39"},
  {name:"C2",X:Q(-1,4),Y:Q(1,4),c:[Q(1),Q(0),Q(1,8),Q(3,8),Q(1,6),Q(5,8),Q(3,8),Q(1,6),Q(0)],cols:[2,3,4,5,6,7],expected:"184320",expectedS:"46/17",expectedExpose:"8/3"}
 ];
 let reports=[];
 for(let s of seeds){
  let p=packets(s.X,s.Y),J=Array.from({length:6},(_,i)=>p.map(x=>x[i]));
  let jets=J.map(r=>r.reduce((a,x,j)=>add(a,mul(x,s.c[j])),zero()));
  if(jets.some(x=>x[0]))throw Error(s.name+" stationary equations");
  let minor=str(det(J.map(r=>s.cols.map(j=>r[j]))));
  if(minor!==s.expected)throw Error(s.name+" minor "+minor);
  let w=mul(Q(4),add(Q(2),add(mul(s.X,s.X),mul(s.Y,s.Y))));
  let S=sub(Q(2),mul(Q(2),div(p[0][0],w))),Se=div(add(mul(w,S),Q(1)),add(w,Q(1,2)));
  if(str(S)!==s.expectedS||str(Se)!==s.expectedExpose)throw Error("CHSH");
  let losing=p.slice(1).map(x=>div(x[0],mul(Q(4),w)));
  if(losing.some(x=>x[0]<=0n))throw Error("losing support");
  reports.push({family:s.name,six_stationary_equations:true,selected_minor:minor,contact_CHSH:str(S),exposing_CHSH:str(Se),all_losing_probabilities_positive:true,losing_probabilities:losing.map(str)});
 }
 const assignments=[{name:"lambda",A:[1,1],B:[1,1],event:0},{name:"n",A:[1,-1],B:[-1,-1],event:1},{name:"beta",A:[-1,1],B:[-1,1],event:7}];
 const events=[[1,1,0,1],[1,-1,0,0],[1,-1,1,1],[-1,1,1,0],[-1,1,0,0],[-1,1,1,1],[-1,-1,0,1],[1,-1,1,0]];
 for(let a of assignments){let probs=events.map(([o,p,x,y])=>Number(a.A[x]===o&&a.B[y]===p)),S=a.A[0]*a.B[0]-a.A[0]*a.B[1]+a.A[1]*a.B[0]+a.A[1]*a.B[1];if(S!==2||probs.some((v,j)=>v!==Number(j===a.event)))throw Error("deterministic neighbor");}
 return{status:"passed",arithmetic:"exact BigInt rationals",seeds:reports,deterministic_neighbors:3,limitation:"Seed calculations corroborate the construction. Universal exposure uses the nonzero contact minors and signs EN1-EN16, as stated in the accompanying proof."};
}
if(typeof process!=="undefined")console.log(JSON.stringify(checkContactSeeds(),null,2));
