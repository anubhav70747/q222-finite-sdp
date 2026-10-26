"use strict";
function checkL12HalfPower(){
 const rows=[
  ["UW",2,0,"W"],["UVX",2,2,"X"],["VZ",0,2,"Z"],["VW",0,2,"W"],
  ["UZ",2,0,"Z"],["UVY",2,2,"Y"],["UVW",2,2,"W"],["Z",0,0,"Z"],
  ["UX",2,0,"X"],["VX",0,2,"X"],["UY",2,0,"Y"],["VY",0,2,"Y"],
  ["X",0,0,"X"],["Y",0,0,"Y"],["sqrt(UV)W",1,1,"W"],["sqrt(UV)Z",1,1,"Z"]
 ];
 if(rows.length!==16||rows.some(x=>x[1]<0||x[1]>2||x[2]<0||x[2]>2))throw Error("outside Z2");
 const reflected=rows.map(([name,i,j,state])=>({name,twice_exponents:[2-j,2-i],state:state==="W"?"Z":state==="Z"?"W":state}));
 if(reflected.some(x=>x.twice_exponents.some(y=>y<0||y>2)))throw Error("inversion outside Z2");
 return{status:"passed",feature_rows:16,feature_convention:"angular exponents are (i/2,j/2)",minimal_grid_containing_displayed_rows:2,rows,reflected,scope:"Checks dictionary and reciprocal-grid inclusion only. Admission and strict-end/compact-uniform positivity are imported from the stated L12 theorem; no arbitrary finite-c or singular-wall claim."};
}
if(typeof process!=="undefined")console.log(JSON.stringify(checkL12HalfPower(),null,2));
