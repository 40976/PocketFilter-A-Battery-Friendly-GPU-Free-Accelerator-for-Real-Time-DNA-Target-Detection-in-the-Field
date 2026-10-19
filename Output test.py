from pynq.overlays.base import BaseOverlay
import time
base = BaseOverlay("base.bit")
# control the LED0
led0 = base.leds[0]
led0.on() # on
time.sleep(1)
led0.off() # off
time.sleep(1)
# repeat
for _ in range(5):
led0.toggle()
time.sleep(0.3)