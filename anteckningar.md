
* U3
_________________________________________________________________________________________________
- d) En tiondels tidssteg på euler jämfört med runge kutta gav fortfarande ett sämre resultat än runge kutta.d

_________________________________________________________________________________________________
- e) 
>> U3

x_slut =

   -1.0607


y_slut =

   -0.4393


felrk4 =

    0.0599


feleuler =

_________________________________________________________________________________________________
f) Rk4 har noggranhetsordningen 4 och euler har noggranhetsordningen 1, enligt förväntan

------------------------
U4)
   Den tar 121 tidssteg. Ploten ser märklig ut (sicksig), men på grund av att felet är så litet misstänker vi att det brusar ur ode45 lösaren, eller låg flyttalsprecision.

------------------------
U5)
   Fall 1: Det blev en spiral eftersom hjulens hastighet ökar med samma acceleration men hjulen har olika hastighet från början vilket får den cirkulera men med växande radie.

   Fall 2: Efterssom accelerationen för hjulen är olika blir spårcirkeln kontinuerligt mindre och slutar till sist i en rand där radien på cirkeln inte minskar mycket men den totala hastigheten ökar.

   Fall 3: Det initiellt långsamma hjulet har en högre acceleration, som leder till att de båda hjulen vid en viss tidpunkt når samma hastighet. Detta följs med att roboten svänger snabbare och snabbare åt andra hållet.

------------------------
U6)
   b) Roboten får radien R = 1.0 för alla t >= 0. Detta beror på att kvoten mellan hastigheterna gör att den rör sig i en perfekt cirkel, även om den kör snabbare och snabbare.
