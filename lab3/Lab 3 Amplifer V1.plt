[FFT of time domain data]
{
   Npanes: 1
   {
      traces: 1 {524290,0,"V(out)"}
      X: ('M',0,10,0,1310710)
      Y[0]: (' ',0,1e-09,20,10)
      Y[1]: (' ',0,-200,40,200)
      Log: 1 2 0
      PltMag: 1
   }
}
[AC Analysis]
{
   Npanes: 1
   {
      traces: 1 {524291,0,"V(out)"}
      X: ('K',0,100,0,10000)
      Y[0]: (' ',3,2.385,0.009,2.475)
      Y[1]: (' ',0,-182,2,-160)
      Units: "ohm" ('M',0,0,0,-1.6e+06,200000,400000)
      Log: 1 0 0
      PltMag: 1
   }
}
[Transient Analysis]
{
   Npanes: 1
   {
      traces: 2 {268959746,0,"V(out)"} {524291,0,"V(in)"}
      X: ('m',0,0,0.01,0.1)
      Y[0]: (' ',1,-2.4,0.4,2.4)
      Y[1]: (' ',0,1e+308,2,-1e+308)
      Volts: (' ',0,0,1,-2.4,0.4,2.4)
      Log: 0 0 0
      PltMag: 1
      PltPhi: 1 0
   }
}
