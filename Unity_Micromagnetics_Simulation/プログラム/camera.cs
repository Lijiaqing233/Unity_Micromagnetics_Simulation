using UnityEngine;
using System.Collections;

//用于控制相机，和UGUI的按钮交互

public class camera : MonoBehaviour
{
    public static bool Kongrotating = true;
    //-----------------  初始参数  -----------------
    //--旋转速度
    private float xSpeed = 125.0f;
    private float ySpeed = 60.0f;

    //中心点物体
    public static Transform Tanks;
    //缩放控制
    private float distance;
    //位置控制
    private Vector3 position;
    //角度控制
    private Quaternion rotation1;
    //相机动画控制：0停止；1旋转指定角度；2旋转到某个角度

    //控制按钮触发的flag
    bool flag_leftR = false;
    bool flag_rightR = false;
    bool flag_upR = false;
    bool flag_downR = false;

    bool flag_leftP = false;
    bool flag_rightP = false;
    bool flag_forwardP = false;
    bool flag_backwardP = false;
    bool flag_upP = false;
    bool flag_downP = false;
    //-----------------  无需改动，中间变量  -----------------
    //旋转角度（二维屏幕）
    public static float x = 90.0f;
    public static float y = 30.0f;
    public static float z = 0.0f;
    //-----------------  公用函数  -----------------
    //相机旋转到某个角度
    //格式化X、Y、Z轴
    public static void FormatXYZ()
    {
        y = y % 360;
        x = x % 360;
        z = z % 360;
    }
    //-----------------  系统函数  -----------------

    

    void Start()
    {
       //初始化中心模块
       Tanks = GameObject.Find("Main Camera").transform;
        //初始化相机位置
        transform.position = new Vector3(-20.8f, 33.8f, 36.9f);
        //初始化相机与中心模型的距离
        //初始化相机角度
        x = 180; y = 45; z = 0;
        //初始化中心物体位置
        Tanks.position = new Vector3(0, 0, 0);
        //初始化中心物体角度
        Tanks.rotation = Quaternion.Euler(0, 0, 0);
        //强制旋转
    }

    void Update()
    {
        if (!Kongrotating) return;
        //如果没有得到中心模块则此脚本无效
        if (Tanks == null) return;
        //初始坐标还原点
        if (Input.GetKey(KeyCode.LeftArrow) | flag_leftR == true) 
        {
            leftR();
        }
        else if (Input.GetKey(KeyCode.RightArrow) | flag_rightR == true)  
        {
            rightR();
        }
        if (Input.GetKey(KeyCode.DownArrow) | flag_downR == true)
        {
            downR();
        }
        else if (Input.GetKey(KeyCode.UpArrow) | flag_upR == true)
        {
            upR();
        }
        if (Input.GetKey(KeyCode.Q) | flag_upP == true)
        {
            UpP();
        }
        else if ( Input.GetKey(KeyCode.E) | flag_downP == true)
        {
            DownP();
        }
        if (Input.GetKey(KeyCode.W) | flag_forwardP == true)
        {
            forwardP();
        }
        else if (Input.GetKey(KeyCode.S) | flag_backwardP == true)
        {
            BackwardP();
        }
        if (Input.GetKey(KeyCode.A) | flag_leftP == true)
        {
            leftP();
        }
        else if (Input.GetKey(KeyCode.D) | flag_rightP == true)
        {
            rightP();
        }
        //格式化 x y z 坐标，使之不可超越-360~360
        FormatXYZ();
        //设置角度
        rotation1 = Quaternion.Euler(y, x, z);
        //设置位置
        position = rotation1 * new Vector3(0.0f, 0.0f, -distance) + Tanks.position;
        //更新角度
        Tanks.rotation = rotation1;
        //更新位置

        //更新角度
        transform.rotation = Tanks.rotation;
    }

    private void upR() {
        y++;
    }

    private void downR() {
        y--;
    }

    private void rightR() {
        x -= 1;
    }

    private void leftR() {
        x += 1;
    }

    private void forwardP()
    {
        transform.Translate(Vector3.forward * 30 * Time.deltaTime);
    }

    private void BackwardP()
    {
        transform.Translate(Vector3.forward * -30 * Time.deltaTime);
    }

    private void rightP()
    {
        transform.Translate(Vector3.right * 60 * Time.deltaTime);
    }

    private void leftP()
    {
        transform.Translate(Vector3.right * -60 * Time.deltaTime);
    }

    private void UpP() {
        transform.Translate(Vector3.up * 30 * Time.deltaTime);
    }

    private void DownP() {
        transform.Translate(Vector3.up * -30 * Time.deltaTime);
    }


    //按键更改UGUI的flag
    public void Enter_leftR() {
        flag_leftR = true;
    }

    public void Exit_leftR() {
        flag_leftR = false;
    }

    public void Enter_rightR()
    {
        flag_rightR = true;
    }

    public void Exit_rightR()
    {
        flag_rightR = false;
    }

    public void Enter_upR()
    {
        flag_upR = true;
    }

    public void Exit_upR()
    {
        flag_upR = false;
    }

    public void Enter_downR()
    {
        flag_downR = true;
    }

    public void Exit_downR()
    {
        flag_downR = false;
    }

    public void Enter_leftP()
    {
        flag_leftP = true;
    }

    public void Exit_leftP()
    {
        flag_leftP = false;
    }

    public void Enter_rightP()
    {
        flag_rightP = true;
    }

    public void Exit_rightP()
    {
        flag_rightP = false;
    }

    public void Enter_forwardP()
    {
        flag_forwardP = true;
    }

    public void Exit_forwardP()
    {
        flag_forwardP = false;
    }

    public void Enter_backwardP()
    {
        flag_backwardP = true;
    }

    public void Exit_backwardP()
    {
        flag_backwardP = false;
    }

    public void Enter_upP()
    {
        flag_upP = true;
    }

    public void Exit_upP()
    {
        flag_upP = false;
    }

    public void Enter_downP()
    {
        flag_downP = true;
    }

    public void Exit_downP()
    {
        flag_downP = false;
    }

}