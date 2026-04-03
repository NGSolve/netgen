#include "array.hpp"
#include "statushandler.hpp"


namespace ngcore
{
  volatile multithreadt multithread;
  
  multithreadt :: multithreadt()
  {
    pause =0;
    testmode = 0;
    redraw = 0;
    drawing = 0;
    terminate = 0;
    running = 0;
    percent = 0;
    task = "";
  }



  
  static Array<std::string> msgstatus_stack(0);
  static Array<double> threadpercent_stack(0);
  static std::string msgstatus = "";


  void ResetStatus()
  {
    SetStatMsg("idle");

    // for (int i = 0; i < msgstatus_stack.Size(); i++)
    // delete msgstatus_stack[i];
    msgstatus_stack.SetSize(0);
    threadpercent_stack.SetSize(0);

    // multithread.task = "";
    multithread.percent = 100.;
  }

  void PushStatus(const std::string& s)
  {
    msgstatus_stack.Append(s);  
    SetStatMsg(s);
    threadpercent_stack.Append(0);
  }
  
  
  void PopStatus()
  {
    if (msgstatus_stack.Size())
      {
        if (msgstatus_stack.Size() > 1)
          // SetStatMsg (*msgstatus_stack.Last());
          SetStatMsg (msgstatus_stack[msgstatus_stack.Size()-2]);
        else
          SetStatMsg ("");
        // delete msgstatus_stack.Last();
        msgstatus_stack.DeleteLast();
        threadpercent_stack.DeleteLast();
        if(threadpercent_stack.Size() > 0)
          multithread.percent = threadpercent_stack.Last();
        else
          multithread.percent = 100.;
      }
    else
      {
        // PrintSysError("PopStatus failed");
        ;
      }
  }
  


  /*
    void SetStatMsgF(const MyStr& s)
    {
    PrintFnStart(s);
    SetStatMsg(s);
    }
  */

  void SetStatMsg(const std::string& s)
  {
    msgstatus = s;
    multithread.task = msgstatus.c_str();  
  }
  
  void SetThreadPercent(double percent)
  {
    multithread.percent = percent;
    if(threadpercent_stack.Size() > 0)
      threadpercent_stack.Last() = percent;
  }


  void GetStatus(std::string & s, double & percentage)
  {
    if(threadpercent_stack.Size() > 0)
      percentage = threadpercent_stack.Last();
    else
      percentage = multithread.percent;
    
    if ( msgstatus_stack.Size() )
      s = msgstatus_stack.Last();
    else
      s = "idle";     
  }
}

