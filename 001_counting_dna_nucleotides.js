var counter = function(instring){
var acount = 0, ccount=0,gcount=0,tcount=0;
for (let index = 0; index < instring.length; index++) {
    if(instring[index]=='A'){
acount++;
    }
    if(instring[index]=='C'){
        ccount++;
    }
    if(instring[index]=='G'){
gcount++;
    }
    if(instring[index]=='T'){
        tcount++;
    }
    
}
    return ""+acount+" "+ccount+" "+gcount+" "+tcount;
}