from kelvin_helmholtz.tools import Kelvin_Helmholtz


def example_of_kelvin_helmholtz():
    print("开始kelvin模块的测试")
    K = Kelvin_Helmholtz()
    print("每步相当于模拟过了时间dt:{}".format(K.dt))
    input("输入回车继续")
    print("print vx")
    print(K.vx)
    print("print vy")
    print(K.vy)
    K.step()
    input("跑了一步. 输入回车继续")
    print("print vx")
    print(K.vx)
    print("print vy")
    print(K.vy)
    input("接下来使用演示功能. 输入回车继续")
    K.display

def main():
    example_of_kelvin_helmholtz()

if __name__ == "__main__":
    main()