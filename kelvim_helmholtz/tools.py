import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec


class Kelvin_Helmholtz(object):
    """"""

    def __init__(self,
                 N=(128, 128),
                 boxsize=(1, 1),
                 courant_fac=0.4,
                 w0=0.1,
                 sigma=0.05/np.sqrt(2),
                 gamma=5/3
                 ):
        self.Nx, self.Ny = N
        self.boxSizeX, self.boxSizeY = boxsize
        self.dx = boxsize[0]/N[0]
        self.dy = boxsize[1]/N[1]
        self.vol = self.dx * self.dy
        self.Y, self.X = np.meshgrid(np.linspace(0.5*self.dy, self.boxSizeY-0.5*self.dy, self.Ny),
                                     np.linspace(0.5*self.dx, self.boxSizeX-0.5*self.dx, self.Nx))
        self.courant_fac = courant_fac
        self.w0 = w0
        self.sigma = sigma
        self.gamma = gamma
        self.rho = 1. + (np.abs(self.Y-0.5) < 0.25)

        # 速度矩阵
        self.vx = -0.5 + (np.abs(self.Y-0.5) < 0.25)
        self.vy = self.w0*np.sin(4*np.pi*self.X) * (np.exp(-(self.Y-0.25)**2/(2 * self.sigma**2)
                                                           ) + np.exp(-(self.Y-0.75)**2/(2*self.sigma**2)))
        self.vz = 0*self.X
        self.P = 0*self.X + 2.5

        # directions for np.roll()
        self.R = -1   # right
        self.L = 1    # left

    def display(self,
                t=(0, 2),
                tOut=0.01,
                useSlopeLimiting=False
                ):
        t, tEnd = t

        def myPlot():
            plt.clf()
            plt.subplot(121)
            plt.imshow(self.rho.T)
            plt.clim(0.8, 2.2)
            ax = plt.gca()
            ax.invert_yaxis()
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            plt.draw()
            plt.subplot(122)
            X = np.linspace(0, 1, 128)
            Y = np.linspace(0, 1, 128)
            plt.streamplot(X, Y, self.vx.T, self.vy.T)
            # plt.quiverkey(q, X=0.3, Y=1.1, U=10,
            #          label='Quiver key, length = 10', labelpos='E')

        # directions for np.roll()
        R = -1   # right
        L = 1    # left

        outputCount = 1

        Mass = self.rho * self.vol
        Momx = self.rho * self.vx * self.vol
        Momy = self.rho * self.vy * self.vol
        Energy = (self.P/(self.gamma-1) + 0.5*self.rho *
                  (self.vx**2+self.vy**2))*self.vol

        while (t < tEnd):
            # get primitive variables
            self.rho = Mass / self.vol
            self.vx = Momx / self.rho / self.vol
            self.vy = Momy / self.rho / self.vol
            self.P = (Energy/self.vol - 0.5*self.rho *
                      (self.vx**2+self.vy**2)) * (self.gamma-1)

            # get time step (CFL)
            dt = self.courant_fac * \
                np.min(np.min([self.dx, self.dy]) /
                       (np.sqrt(self.gamma*self.P/self.rho) + np.sqrt(self.vx**2+self.vy**2)))
            plotThisTurn = False
            if t + dt > outputCount*tOut:
                dt = outputCount*tOut - t
                plotThisTurn = True

            # calculate gradients
            rho_gradx = (np.roll(self.rho, R, axis=0) -
                         np.roll(self.rho, L, axis=0)) / (2.*self.dx)
            rho_grady = (np.roll(self.rho, R, axis=1) -
                         np.roll(self.rho, L, axis=1)) / (2.*self.dy)
            vx_gradx = (np.roll(self.vx, R, axis=0) -
                        np.roll(self.vx, L, axis=0)) / (2.*self.dx)
            vx_grady = (np.roll(self.vx, R, axis=1) -
                        np.roll(self.vx, L, axis=1)) / (2.*self.dy)
            vy_gradx = (np.roll(self.vy, R, axis=0) -
                        np.roll(self.vy, L, axis=0)) / (2.*self.dx)
            vy_grady = (np.roll(self.vy, R, axis=1) -
                        np.roll(self.vy, L, axis=1)) / (2.*self.dy)
            P_gradx = (np.roll(self.P, R, axis=0) -
                       np.roll(self.P, L, axis=0)) / (2.*self.dx)
            P_grady = (np.roll(self.P, R, axis=1) -
                       np.roll(self.P, L, axis=1)) / (2.*self.dy)

            # slope limit gradients
            if useSlopeLimiting:
                rho_gradx = np.maximum(0., np.minimum(
                    1., ((self.rho-np.roll(self.rho, L, axis=0))/self.dx)/(rho_gradx + 1.0e-8*(rho_gradx == 0)))) * rho_gradx
                rho_gradx = np.maximum(0., np.minimum(
                    1., (-(self.rho-np.roll(self.rho, R, axis=0))/self.dx)/(rho_gradx + 1.0e-8*(rho_gradx == 0)))) * rho_gradx
                rho_grady = np.maximum(0., np.minimum(
                    1., ((self.rho-np.roll(self.rho, L, axis=1))/self.dy)/(rho_grady + 1.0e-8*(rho_grady == 0)))) * rho_grady
                rho_grady = np.maximum(0., np.minimum(
                    1., (-(self.rho-np.roll(self.rho, R, axis=1))/self.dy)/(rho_grady + 1.0e-8*(rho_grady == 0)))) * rho_grady
                vx_gradx = np.maximum(0., np.minimum(
                    1., ((self.vx-np.roll(self.vx, L, axis=0))/self.dx) / (vx_gradx + 1.0e-8*(vx_gradx == 0)))) * vx_gradx
                vx_gradx = np.maximum(0., np.minimum(
                    1., (-(self.vx-np.roll(self.vx, R, axis=0))/self.dx) / (vx_gradx + 1.0e-8*(vx_gradx == 0)))) * vx_gradx
                vx_grady = np.maximum(0., np.minimum(
                    1., ((self.vx-np.roll(self.vx, L, axis=1))/self.dy) / (vx_grady + 1.0e-8*(vx_grady == 0)))) * vx_grady
                vx_grady = np.maximum(0., np.minimum(
                    1., (-(self.vx-np.roll(self.vx, R, axis=1))/self.dy) / (vx_grady + 1.0e-8*(vx_grady == 0)))) * vx_grady
                vy_gradx = np.maximum(0., np.minimum(
                    1., ((self.vy-np.roll(self.vy, L, axis=0))/self.dx) / (vy_gradx + 1.0e-8*(vy_gradx == 0)))) * vy_gradx
                vy_gradx = np.maximum(0., np.minimum(
                    1., (-(self.vy-np.roll(self.vy, R, axis=0))/self.dx) / (vy_gradx + 1.0e-8*(vy_gradx == 0)))) * vy_gradx
                vy_grady = np.maximum(0., np.minimum(
                    1., ((self.vy-np.roll(self.vy, L, axis=1))/self.dy) / (vy_grady + 1.0e-8*(vy_grady == 0)))) * vy_grady
                vy_grady = np.maximum(0., np.minimum(
                    1., (-(self.vy-np.roll(self.vy, R, axis=1))/self.dy) / (vy_grady + 1.0e-8*(vy_grady == 0)))) * vy_grady
                P_gradx = np.maximum(0., np.minimum(
                    1., ((self.P-np.roll(self.P, L, axis=0))/self.dx) / (P_gradx + 1.0e-8*(P_gradx == 0)))) * P_gradx
                P_gradx = np.maximum(0., np.minimum(
                    1., (-(self.P-np.roll(self.P, R, axis=0))/self.dx) / (P_gradx + 1.0e-8*(P_gradx == 0)))) * P_gradx
                P_grady = np.maximum(0., np.minimum(
                    1., ((self.P-np.roll(self.P, L, axis=1))/self.dy) / (P_grady + 1.0e-8*(P_grady == 0)))) * P_grady
                P_grady = np.maximum(0., np.minimum(
                    1., (-(self.P-np.roll(self.P, R, axis=1))/self.dy) / (P_grady + 1.0e-8*(P_grady == 0)))) * P_grady

            # extrapolate to cell faces (in time & space)
            rho_prime = self.rho - 0.5*dt * \
                (self.vx * rho_gradx + self.rho * vx_gradx +
                 self.vy * rho_grady + self.rho * vy_grady)
            rho_XL = rho_prime - rho_gradx * self.dx/2.
            rho_XL = np.roll(rho_XL, R, axis=0)
            rho_XR = rho_prime + rho_gradx * self.dx/2.
            rho_YL = rho_prime - rho_grady * self.dy/2.
            rho_YL = np.roll(rho_YL, R, axis=1)
            rho_YR = rho_prime + rho_grady * self.dy/2.
            vx_prime = self.vx - 0.5*dt * (self.vx * vx_gradx + self.vy *
                                           vx_grady + (1/self.rho) * P_gradx)
            vx_XL = vx_prime - vx_gradx * self.dx/2.
            vx_XL = np.roll(vx_XL, R, axis=0)
            vx_XR = vx_prime + vx_gradx * self.dx/2.
            vx_YL = vx_prime - vx_grady * self.dy/2.
            vx_YL = np.roll(vx_YL, R, axis=1)
            vx_YR = vx_prime + vx_grady * self.dy/2.
            vy_prime = self.vy - 0.5*dt * (self.vx * vy_gradx + self.vy *
                                           vy_grady + (1/self.rho) * P_grady)
            vy_XL = vy_prime - vy_gradx * self.dx/2.
            vy_XL = np.roll(vy_XL, R, axis=0)
            vy_XR = vy_prime + vy_gradx * self.dx/2.
            vy_YL = vy_prime - vy_grady * self.dy/2.
            vy_YL = np.roll(vy_YL, R, axis=1)
            vy_YR = vy_prime + vy_grady * self.dy/2.
            P_prime = self.P - 0.5*dt * \
                (self.gamma*self.P * (vx_gradx + vy_grady) +
                 self.vx * P_gradx + self.vy * P_grady)
            P_XL = P_prime - P_gradx * self.dx/2.
            P_XL = np.roll(P_XL, R, axis=0)
            P_XR = P_prime + P_gradx * self.dx/2.
            P_YL = P_prime - P_grady * self.dy/2.
            P_YL = np.roll(P_YL, R, axis=1)
            P_YR = P_prime + P_grady * self.dy/2.

            # compute star (averaged) states
            rho_Xstar = 0.5*(rho_XL + rho_XR)
            rho_Ystar = 0.5*(rho_YL + rho_YR)
            momx_Xstar = 0.5*(rho_XL * vx_XL + rho_XR * vx_XR)
            momx_Ystar = 0.5*(rho_YL * vx_YL + rho_YR * vx_YR)
            momy_Xstar = 0.5*(rho_XL * vy_XL + rho_XR * vy_XR)
            momy_Ystar = 0.5*(rho_YL * vy_YL + rho_YR * vy_YR)
            en_Xstar = 0.5*(P_XL/(self.gamma-1)+0.5*rho_XL * (vx_XL**2+vy_XL **
                                                              2) + P_XR/(self.gamma-1)+0.5*rho_XR * (vx_XR**2+vy_XR**2))
            en_Ystar = 0.5*(P_YL/(self.gamma-1)+0.5*rho_YL * (vx_YL**2+vy_YL **
                                                              2) + P_YR/(self.gamma-1)+0.5*rho_YR * (vx_YR**2+vy_YR**2))

            P_Xstar = (self.gamma-1)*(en_Xstar-0.5 *
                                      (momx_Xstar**2+momy_Xstar**2)/rho_Xstar)
            P_Ystar = (self.gamma-1)*(en_Ystar-0.5 *
                                      (momx_Ystar**2+momy_Ystar**2)/rho_Ystar)

            # compute fluxes (local Lax-Friedrichs/Rusanov)
            flux_rho_X = momx_Xstar
            flux_rho_Y = momy_Ystar
            flux_momx_X = momx_Xstar**2/rho_Xstar + P_Xstar
            flux_momx_Y = momy_Ystar * momx_Ystar/rho_Ystar
            flux_momy_X = momx_Xstar * momy_Xstar/rho_Xstar
            flux_momy_Y = momy_Ystar**2/rho_Ystar + P_Ystar
            flux_en_X = (en_Xstar+P_Xstar) * momx_Xstar/rho_Xstar
            flux_en_Y = (en_Ystar+P_Ystar) * momy_Ystar/rho_Ystar

            C = np.sqrt(self.gamma*P_XL/rho_XL) + np.abs(vx_XL)
            C = np.maximum(C, np.sqrt(self.gamma*P_XR/rho_XR) + np.abs(vx_XR))
            C = np.maximum(C, np.sqrt(self.gamma*P_YL/rho_YL) + np.abs(vy_YL))
            C = np.maximum(C, np.sqrt(self.gamma*P_YR/rho_YR) + np.abs(vy_YR))

            flux_rho_X = flux_rho_X - C * 0.5 * (rho_XL - rho_XR)
            flux_rho_Y = flux_rho_Y - C * 0.5 * (rho_YL - rho_YR)
            flux_momx_X = flux_momx_X - C * 0.5 * \
                (rho_XL * vx_XL - rho_XR * vx_XR)
            flux_momx_Y = flux_momx_Y - C * 0.5 * \
                (rho_YL * vx_YL - rho_YR * vx_YR)
            flux_momy_X = flux_momy_X - C * 0.5 * \
                (rho_XL * vy_XL - rho_XR * vy_XR)
            flux_momy_Y = flux_momy_Y - C * 0.5 * \
                (rho_YL * vy_YL - rho_YR * vy_YR)
            flux_en_X = flux_en_X - C * 0.5 * (P_XL/(self.gamma-1)+0.5*rho_XL * (
                vx_XL**2+vy_XL**2) - (P_XR/(self.gamma-1)+0.5*rho_XR * (vx_XR**2+vy_XR**2)))
            flux_en_Y = flux_en_Y - C * 0.5 * (P_YL/(self.gamma-1)+0.5*rho_YL * (
                vx_YL**2+vy_YL**2) - (P_YR/(self.gamma-1)+0.5*rho_YR * (vx_YR**2+vy_YR**2)))

            # update solution
            Mass = Mass - dt * self.dy * flux_rho_X
            Mass = Mass + dt * self.dy * np.roll(flux_rho_X, L, axis=0)
            Mass = Mass - dt * self.dx * flux_rho_Y
            Mass = Mass + dt * self.dx * np.roll(flux_rho_Y, L, axis=1)
            Momx = Momx - dt * self.dy * flux_momx_X
            Momx = Momx + dt * self.dy * np.roll(flux_momx_X, L, axis=0)
            Momx = Momx - dt * self.dx * flux_momx_Y
            Momx = Momx + dt * self.dx * np.roll(flux_momx_Y, L, axis=1)
            Momy = Momy - dt * self.dy * flux_momy_X
            Momy = Momy + dt * self.dy * np.roll(flux_momy_X, L, axis=0)
            Momy = Momy - dt * self.dx * flux_momy_Y
            Momy = Momy + dt * self.dx * np.roll(flux_momy_Y, L, axis=1)
            Energy = Energy - dt * self.dy * flux_en_X
            Energy = Energy + dt * self.dy * np.roll(flux_en_X, L, axis=0)
            Energy = Energy - dt * self.dx * flux_en_Y
            Energy = Energy + dt * self.dx * np.roll(flux_en_Y, L, axis=1)

            # advance time
            t += dt

            # plot the solution at regular time intervals
            if plotThisTurn:
                print(t)
                myPlot()
                plt.pause(0.001)
                outputCount += 1

        plt.show()


def main():
    K = Kelvin_Helmholtz()
    w0 = 0.1
    sigma = 0.05/np.sqrt(2)
    rho = 1. + (np.abs(K.Y-0.5) < 0.25)
    vx = -0.5 + (np.abs(K.Y-0.5) < 0.25)
    vy = w0*np.sin(4*np.pi*K.X) * (np.exp(-(K.Y-0.25)**2/(2 * sigma**2)
                                          ) + np.exp(-(K.Y-0.75)**2/(2*sigma**2)))
    vz = 0*K.X
    P = 0*K.X + 2.5

    K.display()


if __name__ == "__main__":
    main()
