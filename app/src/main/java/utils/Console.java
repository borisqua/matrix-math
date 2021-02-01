package utils;

import java.io.IOException;
import java.util.concurrent.TimeUnit;

public class Console {
  
  public static void clear() {
    try {
      String os = System.getProperty("os.name");
    if (os.contains("Windows")) {
        new ProcessBuilder("cmd", "/c", "cls").inheritIO().start().waitFor();
    } else if (os.contains("Linux")) {
      System.out.print("\033[H\033[2J");
    } else {
      Runtime.getRuntime().exec("clear");
    }
    } catch (InterruptedException | IOException ignored) {
    }
  }
  
  public static void waitAndThenClean(int fps) {
    try {
      TimeUnit.MILLISECONDS.sleep(1000 / fps);
      clear();
    } catch (InterruptedException ignored) {
    }
  }
  
  public static void showProgress(int learningIteration, int iterationsLimit, int step, String caption) {
    if (learningIteration % step != 0) {
      return;
    }
    String unicodeBlockSymbol = (new String(new char[]{'\u2588'}));
    clear();
    double percentsDone = 100D * learningIteration / iterationsLimit;
    System.out.printf(caption + " %2.2f%% done\n", percentsDone);
    System.out.println(unicodeBlockSymbol.repeat((int) (50 * percentsDone / 100)));
  }
  
}
