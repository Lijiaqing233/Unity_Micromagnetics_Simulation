using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;




public class UiBase : MonoBehaviour
{
  

    public InputField Nx;
    public InputField Ny;
    public InputField Nz;
    public InputField Hxc;
    public InputField Hyc;
    public InputField Hzc;
    public InputField Dx;
    public InputField Dy;
    public InputField Dz;
    public InputField Dt;
    public InputField Al;
    public void FindUI()
    {
        Nx = GameObject.Find("Nx").GetComponent<InputField>();
        Ny = GameObject.Find("Ny").GetComponent<InputField>();
        Nz = GameObject.Find("Nz").GetComponent<InputField>();
        Hxc = GameObject.Find("Hxc").GetComponent<InputField>();
        Hyc = GameObject.Find("Hyc").GetComponent<InputField>();
        Hzc = GameObject.Find("Hzc").GetComponent<InputField>();
        Dt = GameObject.Find("Dt").GetComponent<InputField>();
        Al = GameObject.Find("Alpha").GetComponent<InputField>();

    }

    public void UpdateDisplayText(Parameter[] parameter)
    {

    }

    




}
