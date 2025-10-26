module m_derivs
    use, intrinsic :: iso_fortran_env, only: i32 => int32, r64 => real64
    use m_hydraulics, only: riverhydraulics_type
    use m_integration_data, only: &
        saveheatfluxtribs, saveheatfluxadvecdisp, saveheatfluxjsnt, saveheatfluxlongat, saveheatfluxback, &
        saveheatfluxconv, saveheatfluxevap, saveheatfluxjsed, saveheatfluxjhyporheic, &
        os, phitsave, savedofluxheadwater, saveco2fluxheadwater, savedofluxtribs, savedofluxadvecdisp, &
        saveco2fluxtribs, saveco2fluxadvecdisp, phinsave, phipsave, phicsave, philsave, phitotalsave , &
        savebotalgphoto, savebotalgresp, savebotalgdeath, savebotalgnetgrowth, &
        sodpr, jnh4pr, jno3pr, jch4pr, jsrppr, csodpr, &
        diagfluxdo, diagfluxcbod, diagfluxnh4, diagfluxno3, diagfluxsrp, diagfluxic, savedofluxreaer, &
        savedofluxcbodfast, savedofluxcbodslow, savedofluxnitrif,  savedofluxphytoresp, &
        savedofluxphytophoto, savedofluxbotalgresp, savedofluxbotalgphoto, savedofluxsod, savedofluxcod, &
        saveco2fluxreaer, saveco2fluxcbodfast, saveco2fluxcbodslow, saveco2fluxphytoresp, saveco2fluxphytophoto, &
        saveco2fluxbotalgresp, saveco2fluxbotalgphoto, saveco2fluxsod, hypofluxdo, hypofluxcbod, hypofluxnh4, &
        hypofluxno3, hypofluxsrp, hypofluxic
    use m_light_heat, only: lightheat, heatbalance, lightextinction
    use m_phsolve, only: chemrates, modfp2, ph_solver
    use m_solar_calc, only: solar_type, solarcalc
    use m_sourcein, only: load, sourcescalc
    use m_system_params, only: system_params_t
    use m_downstream_boundary, only: downstream_boundary_t, instanteousdownstreamboundary
    use m_meteorology, only: meteorology_t, instanteousmeteo
    use m_output, only: outdata_t, output
    use m_oxygen, only: oxygen_inhibition_and_enhancement, oxygen_saturation
    use m_rates, only: rates_t
    use m_tempadjust, only: temp_adjust
    use m_upstream_boundary, only: upstream_boundary_t, instanteousheadwater
    use m_constants, only: nv, nl, adam, bdam, cpw, rhow
    use wq_idx
    use m_sedcalcnew, only: sedcalcnumnew
    implicit none
    private
    public :: derivs

contains

    subroutine derivs(nr, meteo, solar, hw, db, hydrau, sys, te, c, inb, ipb, &
        rates, dte, dc, dinb, dipb, t)

        implicit none
        integer(i32), intent(in) :: nr
        type(meteorology_t) :: meteo
        type(solar_type) solar !solar radiation
        type(upstream_boundary_t) hw !headwater
        type(downstream_boundary_t) db !downstream boundary
        type(riverhydraulics_type) hydrau !channel dimensions, hydraulics, physical characters
        type(system_params_t) sys
        real(r64), dimension(:), pointer :: inb, ipb
        !real(r64), dimension(:), pointer :: phitotalsave, phitsave, philsave, phinsave, phipsave, phicsave !gp 20-oct-04
        !gp 27-oct-04 real(r64), dimension(:,:), pointer :: te, c
        real(r64), dimension(:,:), pointer :: te !gp 27-oct-04
        real(r64), dimension(:,:,:), pointer :: c !gp add dimension for nl
        type(rates_t) rates
        !gp 27-oct-04 real(r64), intent(out):: dte(nr, nl), dc(nr, nv), dinb(nr), dipb(nr)
        real(r64), intent(out):: dte(nr, nl), dc(nr, nv, nl), dinb(nr), dipb(nr) !gp
        real(r64), intent(in):: t !time

        integer(i32) i, k
        real(r64) :: evap, back, flux
        real(r64) longat, conv
        !temperature corrected rates

        !gp 03-apr-08
        !real(r64) khct(nr), kdcst(nr), kdct(nr), kgaft(nr), kdeaft(nr), kreaft(nr), kexaft(nr)
        real(r64) khct(nr), kdcst(nr), kdct(nr), kgaft(nr), kdeaft(nr), krea1ft(nr), krea2ft(nr), kexaft(nr)

        real(r64) khnt(nr), khpt(nr), knt(nr), kdtt(nr), kpatht(nr)
        real(r64) kit(nr), kat(nr), kgat(nr), kdeat(nr), kreat(nr), vdit(nr), kact(nr)

        !gp 30-nov-04
        real(r64) kgent(nr) !gp 30-nov-04 generic constituent decay

        real(r64) :: ph=0
        real(r64) fcarb, fnitr, fdenitr, frespp, frespb
        real(r64) ke, phip, phint
        real(r64) phil, num, den
        real(r64) alpha0, alpha1

        real(r64) prefam, prefamf
        real(r64) sod, jch4, jnh4, jno3, jpo4

        real(r64) jcin, jnin, jpin
        real(r64) ow, nh3w, no3w, po4w, tw
        real(r64) jamm, jnitr, jmeth, jmethg, jphos

        real(r64) botlight
        real(r64) dam, dropft, defa, defb, iat(nr) !changed by hua

        !ph
        real(r64) k1, k2, kw, kh, co2sat

        !gp 01-nov-07
        !real(r64) hh, alp0, alp1, chco3co3
        real(r64) hh, alp0, alp1, alp2, chco3co3

        real(r64) :: csod =0

        real(r64) detrdiss, detrsettl
        real(r64) cbodshydr, cbodsoxid, cbodfoxid
        real(r64) orgnhydr, nh4nitrif, denitr
        real(r64) orgphydr
        real(r64) orgnsettl, orgpsettl, inorgpsettl
        real(r64) inorgsettl
        real(r64) oxreaer
        real(r64) phytophoto, phytoresp, phytodeath, phytosettl
        real(r64) botalgphoto, botalgresp, botalgdeath, botalgexc
        real(r64) botalguptaken, botalguptakep
        real(r64) :: ninb = 0, nipb =0, fbnb =0, fbpb =0
        real(r64) :: phic =0
        real(r64) :: plim, nlim

        !gp 15-nov-04 integer(i32) :: hco3usef

        !pathogens
        real(r64) :: ksol=0

        !gp 27-oct-04 real(r64) phs(0:nr), k1s(0:nr), k2s(0:nr), khs(0:nr) !new 08/29/04 for less ph solver be called
        real(r64) phs(0:nr, nl), k1s(0:nr, nl), k2s(0:nr, nl), khs(0:nr, nl) !gp add dim for nl

        !gp 15-nov-04 move ehyporheic to hydrau
        !gp 20-oct-04 sediment and hyporheic heat flux variables
        !gp real(r64) jsed(nr), jhyporheic(nr), ehyporheic(nr)
        real(r64) jsed(nr), jhyporheic(nr)

        !gp 03-nov-04 hyporheic kinetics variables
        real(r64) kgaht(nr), fcarbh, heterogrow

        !gp 15-nov-04 additional variables for level 2 hyporheic biofilm kinetics
        real(r64) kreaht(nr), kdeaht(nr), heteroresp, heterodeath, prefamh

        !gp 08-dec-04 cod oxidation if generic constituent is used as cod
        real(r64) codoxid

        !gp 03-dec-09 variables for fraction of ionized nh4+ and phosphate speciation
        real(r64) fi, kamm
        real(r64) kpo41, kpo42,kpo43
        real(r64) dpo
        real(r64) fpo41, fpo42,fpo43
        real(r64) k1nh3, kgnh3, henh3, vnh3, naus, nh3gas
        real(r64) pcharge, uw10

        dte =0; dc =0; dinb =0; dipb =0
        !gp move kawind to this sub to use interpolated hourly wind speed
        !real(r64) kawind

        sys%tday = t - int(t)

        call instanteousmeteo(nr, t, meteo)

        !scc add the following loop to calculate ph and save results 08/29/04
        do i=0, nr

            call ph_solver(sys%imethph, ph, c(i, IDX_TOTAL_INORGANIC_C, 1), te(i, 1), &
                c(i, IDX_ALKALINITY, 1), c(i, IDX_CONDUCTIVITY, 1))

            phs(i, 1) = ph
            phs(i, 2) = 0 !not used unless hyporheic wq is simulated
            call chemrates(te(i, 1), k1, k2, kw, kh, c(i, IDX_CONDUCTIVITY, 1))
            k1s(i, 1) = k1; k2s(i, 1) = k2; khs(i, 1) = kh
            k1s(i, 2) = 0; k2s(i, 2) = 0; khs(i, 2) = 0 !gp 02-nov-04 end of new block of code
        end do

        ! uw(0) = uw(1)

        !gp evaluate point source sine functions and distribute loads to reaches for the current time t
        call sourcescalc(t, nr)

        call solarcalc(nr, solar, meteo, hydrau, sys) !solarcalc(tday, i)

        !gp 27-oct-04 call instanteousheadwater(hw, t, te(0,1), c(0,:), ph)
        !gp call instanteousdownstreamboundary(db, t, te(nr,1), c(nr,:), te(nr+1,1), c(nr+1,:))
        call instanteousheadwater(hw, t, te(0,1), c(0,:,1), ph)

        !gp 12-jan-06
        !write(10,*) 'done thru call instanteousheadwater'

        call instanteousdownstreamboundary(db, t, te(nr,1), c(nr,:,1), te(nr+1,1), c(nr+1,:,1)) !gp 27-oct-04 end new block

        !gp 12-jan-06
        !write(10,*) 'done thru call instanteousdownstreamboundary'
        !write(10,'(32f13.4)') te(nr,1), (c(nr, k, 1), k=1, nv-2), te(nr+1,1), (c(nr+1, k, 1), k=1, nv-2)


        !'
        !' ------------------------------------------------------------------------------------
        !' --- heat transport and flux derivatives for the water column and sediment layers ---
        !' ------------------------------------------------------------------------------------
        !'

        !heat transport derivatives

        do i = 1, nr
            ! dte(i, 1) = hydrau%reach(i-1)%qcmd * te(i - 1, 1) * rhow * cpw &
            ! - hydrau%reach(i)%qcmd * te(i, 1) * rhow * cpw &
            ! + hydrau%reach(i)%qptcmd * load(i)%te * rhow * cpw &
            ! - hydrau%reach(i)%qptacmd * te(i, 1) * rhow * cpw &
            ! + hydrau%reach(i-1)%epcmd * (te(i - 1, 1) - te(i, 1)) * rhow * cpw &
            ! + hydrau%reach(i)%epcmd * (te(i + 1, 1) - te(i, 1)) * rhow * cpw
            dte(i,1) = hydrau%reach(i-1)%qcmd * te(i - 1, 1) * rhow * cpw
            dte(i,1) = dte(i,1) - hydrau%reach(i)%qcmd * te(i, 1) * rhow * cpw
            dte(i,1) = dte(i,1) + hydrau%reach(i)%qptcmd * load(i)%te * rhow * cpw
            dte(i,1) = dte(i,1) - hydrau%reach(i)%qptacmd * te(i, 1) * rhow * cpw
            dte(i,1) = dte(i,1) + hydrau%reach(i-1)%epcmd * (te(i - 1, 1) - te(i, 1)) * rhow * cpw
            dte(i,1) = dte(i,1) + hydrau%reach(i)%epcmd * (te(i + 1, 1) - te(i, 1)) * rhow * cpw

            !'gp 05-jul-05 (output fluxes in cal/cm^2/d)
            saveheatfluxtribs(i) = ((hydrau%reach(i)%qptcmd * load(i)%te &
                - hydrau%reach(i)%qptacmd * te(i, 1)) &
                * rhow * cpw * 100 / hydrau%reach(i)%ast)
            saveheatfluxadvecdisp(i) = ((hydrau%reach(i-1)%qcmd * te(i - 1, 1) &
                - hydrau%reach(i)%qcmd * te(i, 1) &
                + hydrau%reach(i-1)%epcmd * (te(i - 1, 1) - te(i, 1)) &
                + hydrau%reach(i)%epcmd * (te(i + 1, 1) - te(i, 1))) &
                * rhow * cpw * 100 / hydrau%reach(i)%ast)

        end do

        do i = 1, nr

            !gp 16-jul-08
            !call heatbalance(evap, conv, back, longat, lightheat, te(i,1), meteo%ta(i), &
            ! meteo%td(i), meteo%cc(i), meteo%uw(i), hydrau%reach(i)%ast)
            call heatbalance(evap, conv, back, longat, lightheat, te(i,1), meteo%ta(i), &
                meteo%td(i), meteo%cc(i), meteo%uw(i), hydrau%reach(i)%ast, hydrau%reach(i)%skop)

            flux = solar%jsnt(i) + longat - back - conv - evap !gp changed jsnt to jsnt(i)
            dte(i, 1) = dte(i, 1) + flux * hydrau%reach(i)%ast / 100.0_r64
            ! if (i==1) write(8,*) meteo%ta(i),meteo%td(i) !, solar%jsnt(i), longat, back, conv, evap

            !'gp 01-jul-05 save heat flux terms (cal/cm2/d) for each reach for this time step for output
            saveheatfluxjsnt(i) = solar%jsnt(i)
            saveheatfluxlongat(i) = longat
            saveheatfluxback(i) = -back
            saveheatfluxconv(i) = -conv
            saveheatfluxevap(i) = -evap

        end do

        do i = 1 , nr
            !'gp new sediment-water heat flux from sediments into the water in units of cal/cm2/day
            jsed(i) = (te(i, 2) - te(i, 1)) * 86400 * 2 * hydrau%reach(i)%sedthermcond / &
                hydrau%reach(i)%hsedcm
            !units of (cal/cm2/day), note that sedthermcond is in units of (cal/sec) per (cm deg c)

            !'gp 05-jul-05 save heat flux terms (cal/cm2/d) for each reach for this time step for output
            saveheatfluxjsed(i) = jsed(i)

            dte(i, 1) = dte(i, 1) + jsed(i) * hydrau%reach(i)%ast/ 100 !'units of (cal/day) * (m3/cm3)
            dte(i, 2) = -jsed(i) * hydrau%reach(i)%ast / 100 !'units of (cal/day) * (m3/cm3)
            select case (sys%simhyporheicwq)
              case ('Level 1', 'Level 2') !gp 15-nov-04
                !'gp new hyporheic exchange flux in units of (cal/day) * (m3/cm3)
                jhyporheic(i) = hydrau%reach(i)%ehyporheiccmd * &
                    (te(i, 2) - te(i, 1)) * rhow * cpw * 100 / hydrau%reach(i)%ast !'units of (cal/cm2/day)

                !'gp 05-jul-05 save heat flux terms (cal/cm2/d) for each reach for this time step for output
                saveheatfluxjhyporheic(i) = jhyporheic(i)

                dte(i, 1) = dte(i, 1) + jhyporheic(i) * hydrau%reach(i)%ast / 100 !'units of (cal/day) * (m3/cm3)
                dte(i, 2) = dte(i, 2) - jhyporheic(i) * hydrau%reach(i)%ast / 100 !'units of (cal/day) * (m3/cm3)
            end select
        end do
        do i = 1, nr
            !'gp units for dte on lhs in the next original line of code is deg c per day
            dte(i, 1) = dte(i, 1) / &
                hydrau%reach(i)%vol / rhow / cpw !'units of deg c per day
            !'gp sediment-water heat flux term dts in deg c per day
            dte(i, 2) = dte(i, 2) / &
                (hydrau%reach(i)%ast * hydrau%reach(i)%hsedcm / 100) / &
                hydrau%reach(i)%rhocpsed !'units of deg c per day
        end do


        !gp 03-feb-05 start of bypass derivs for water quality variables
        !unless 'all' state variables are being simulated
        !
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        !gp 29-dec-09
        !if (sys%statevariables == "all") then
        if (sys%statevariables /= "Temperature") then

            !x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
            !x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x



            !gp 02-nov-04 os(0) = oxsat(te(0, 1), hydrau%reach(0)%elev)
            os(0, 1) = oxygen_saturation(te(0, 1), hydrau%reach(0)%elev)
            select case (sys%simhyporheicwq)
              case ('Level 1', 'Level 2') !gp 15-nov-04
                os(0, 2) = oxygen_saturation(te(0, 2), hydrau%reach(0)%elev)
              case default
                os(0, 2) = 0 !not used in calculation
            end select !gp 02-nov-04 end new block

            !'
            !' --------------------------------------------------------------------------
            !' --- temperature adjustment of rates and constants for the water column ---
            !' --------------------------------------------------------------------------
            !'

            !gp 03-apr-08
            !call tempadjust(rates, hydrau,meteo, nr, te, khct, kdcst, kdct, kgaft, kdeaft, &
            ! kreaft, kexaft, khnt, khpt, knt, kdtt, kpatht, kit, kat, kgat, kdeat, &
            ! kreat, vdit, kact, kgent) !gp 30-nov-04 add kgent
            call temp_adjust(rates, hydrau,meteo, nr, te, khct, kdcst, kdct, kgaft, kdeaft, &
                krea1ft, krea2ft, kexaft, khnt, khpt, knt, kdtt, kpatht, kit, kat, kgat, kdeat, &
                kreat, vdit, kact, kgent) !gp 30-nov-04 add kgent


            !'gp 20-oct-04 temp limitation factor for output
            do i = 1 , nr
                phitsave(i) = rates%tkgaf ** (te(i, 1) - 20)
            end do

            do i=1 ,nr
                !gp 27-oct-04 os(i)=oxsat(te(i,1), hydrau%reach(i)%elev)
                os(i, 1)=oxygen_saturation(te(i, 1), hydrau%reach(i)%elev) !gp 27-oct-04
            end do

            !'
            !' ---------------------------------------------------------------------------
            !' ---------------------------------------------------------------------------
            !' --- ---
            !' --- water column water quality deriviatives ---
            !' --- note that units are mass/time in the main i loop ---
            !' --- and at the the derivs are divided by vol(i) to get mass/volume/time ---
            !' --- or area to get mass/area/time ---
            !' --- ---
            !' ---------------------------------------------------------------------------
            !' ---------------------------------------------------------------------------
            !'

            !transport and small dams
            do i = 1, nr
                do k = 1, nv - 1
                    select case (k)

                        !gp 31-mar-05 ct is variable nv-1 (debug alkalinity deriv)
                      case (3,nv-1)

                        dropft = hydrau%reach(i -1)%drop * 3.281_r64 !scc
                        dam = 1.0_r64 / (1 + 0.116_r64 * adam * bdam * dropft * (1 - 0.034_r64 * dropft) * &
                            (1.0_r64 + 0.046_r64 * te(i - 1, 1)))
                        if (k == 3) then
                            defa = os(i - 1, 1) - c(i - 1, IDX_DISSOLVED_OXYGEN, 1)
                            defb = dam * defa
                            dc(i, k, 1) = hydrau%reach(i-1)%qcmd * (os(i-1, 1) - defb) !gp 27-oct-04 end new block

                            !'gp 05-jul-05 (save do fluxes in go2/m^2/d)
                            savedofluxheadwater(i) = (hydrau%reach(i-1)%qcmd * (os(i - 1, 1) - defb) &
                                - hydrau%reach(i-1)%qcmd * c(i - 1, k, 1)) / hydrau%reach(i)%ast

                            !gp 31-mar-05 ct is variable nv-2
                            !elseif (k == 15) then
                        elseif (k == nv-1) then

                            co2sat = khs(i-1, 1) * rates%pco2
                            hh = 10.0_r64 ** (-phs(i-1, 1))
                            alp0 = hh * hh / (hh ** 2.0 + k1s(i-1, 1) * hh + k1s(i-1, 1) * k2s(i-1, 1))
                            chco3co3 = (1.0_r64 - alp0) * c(i-1, k, 1)
                            defa = co2sat - alp0 * c(i-1, k, 1)
                            defb = dam * defa
                            dc(i, k, 1) = hydrau%reach(i-1)%qcmd * (chco3co3 + co2sat - defb) !gp 27-oct-04 end new block

                            !'gp 05-jul-05 (save co2 fluxes in gc/m^2/d)
                            saveco2fluxheadwater(i) = (hydrau%reach(i-1)%qcmd * (chco3co3 + co2sat - defb) &
                                - hydrau%reach(i-1)%qcmd * c(i - 1, k, 1)) &
                                / hydrau%reach(i)%ast / rates%rccc

                        end if
                      case default !else
                        dc(i, k, 1) = hydrau%reach(i-1)%qcmd * c(i - 1, k, 1) !gp 27-oct-04 end new block
                    end select
                    dc(i, k, 1) = dc(i, k, 1) - hydrau%reach(i)%qcmd * c(i, k, 1) &
                        + hydrau%reach(i)%qptcmd * load(i)%c(k) &
                        - hydrau%reach(i)%qptacmd * c(i, k, 1) &
                        + hydrau%reach(i-1)%epcmd * (c(i - 1, k, 1) - c(i, k, 1)) &
                        + hydrau%reach(i)%epcmd * (c(i + 1, k, 1) - c(i, k, 1)) !gp 27-oct-04 end new block

                    !'gp 05-jul-05 (save net trib/gw and advec/disp do fluxes in go2/m^2/d and co2 fluxes in gc/m^2/d)
                    select case (k)
                      case (3)
                        savedofluxtribs(i) = (hydrau%reach(i)%qptcmd * load(i)%c(k) &
                            - hydrau%reach(i)%qptacmd * c(i, k, 1)) &
                            / hydrau%reach(i)%ast
                        savedofluxadvecdisp(i) = savedofluxheadwater(i) &
                            + (-hydrau%reach(i)%qcmd * c(i, k, 1) &
                            + hydrau%reach(i-1)%epcmd * (c(i - 1, k, 1) - c(i, k, 1)) &
                            + hydrau%reach(i)%epcmd * (c(i + 1, k, 1) - c(i, k, 1))) &
                            / hydrau%reach(i)%ast
                      case (nv - 1)
                        saveco2fluxtribs(i) = (hydrau%reach(i)%qptcmd * load(i)%c(k) &
                            - hydrau%reach(i)%qptacmd * c(i, k, 1)) &
                            / hydrau%reach(i)%ast / rates%rccc
                        saveco2fluxadvecdisp(i) = saveco2fluxheadwater(i) &
                            + (-hydrau%reach(i)%qcmd * c(i, k, 1) &
                            + hydrau%reach(i-1)%epcmd * (c(i - 1, k, 1) - c(i, k, 1)) &
                            + hydrau%reach(i)%epcmd * (c(i + 1, k, 1) - c(i, k, 1))) &
                            / hydrau%reach(i)%ast / rates%rccc
                    end select

                end do
            end do

            ! ----------------
            ! --- kinetics ---
            ! ----------------

            do i = 1, nr

                ! --- ph dependent variables for later use in deriv calcs in this i loop

                hh = 10.0_r64 ** (-phs(i, 1))
                alp0 = hh * hh / (hh ** 2.0_r64 + k1s(i, 1) * hh + k1s(i, 1) * k2s(i, 1))
                alp1 = k1s(i, 1) * hh / (hh ** 2.0_r64 + k1s(i, 1) * hh + k1s(i, 1) * k2s(i, 1)) !fraction of ct as hco3-
                alp2 = k1s(i, 1) * k2s(i, 1) / (hh ** 2.0_r64 + k1s(i, 1) * hh + k1s(i, 1) * k2s(i, 1)) !fraction of ct as co3--

                !'gp 03-dec-09
                !'fraction of total ammonia that is ionized ammonia nh4+ (fi)
                kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (te(i, 1) + 273.15_r64)))
                !equilibrium coeff for ammonia dissociation nh4+ = nh3(aq) + h+
                fi = hh / (hh + kamm) !'fraction of ionized nh4+ = [nh4+ / (nh4+ + nh3)]
                !'fraction of srp that is h2po4- (fpo41), hpo4-- (fpo42), and po4--- (fpo43)
                kpo41 = 10.0_r64 ** (-2.15_r64) !'= [h+][h2po4-]/[h3po4] for eqn h3po4 = h+ + h2po4-
                kpo42 = 10.0_r64 ** (-7.2_r64) !'= [h+][hpo4--]/[h2po4-] for eqn h2po4- = h+ + hpo4--
                kpo43 = 10.0_r64 ** (-12.35_r64) !'= [h+][po4---]/[hpo4--] for eqn hpo4-- = h+ + po4---
                dpo = 1.0_r64 / (hh ** 3.0_r64 + kpo41 * hh ** 2.0_r64 + kpo41 * kpo42 * hh + kpo41 * kpo42 * kpo43)
                !intermediate calc of denominator
                fpo41 = kpo41 * hh ** 2.0_r64 * dpo !'fraction of phosphate as h2po4-
                fpo42 = kpo41 * kpo42 * hh * dpo !'fraction of phosphate as hpo4--
                fpo43 = kpo41 * kpo42 * kpo43 * dpo !'fraction of phosphate as po4---
                !'stoichiometric conversion factors for alkalinity (some are ph dependent) (ralkbn & ralkbp calc in classstoch.f90)
                !ralkbn = 1# / 14.0067 / 1000# / 1000# !'eqh+/l per mgn/m^3 (1eqh+/molen / 14gn/molen / 1000mgn/gn / 1000l/m^3)
                !ralkbp = (fpo41 + 2# * fpo42 + 3# * fpo43) / 30.973762 / 1000# / 1000# !'eqh+/l per mgp/m^3 (eqh+/molep / 31gp/molep / 1000mgp/gp / 1000l/m^3)
                pcharge = (fpo41 + 2.0_r64 * fpo42 + 3.0_r64 * fpo43) !'average charge of p species to multiply by ralkbp
                !'ammonia gas transfer at air/water interface (eqns from chapra et al qual2k manual ver 2.11b8)
                k1nh3 = 1.171_r64 * kat(i) * hydrau%reach(i)%depth !'liquid film exchange coefficient for nh3 (m/d)
                uw10 = (10.0_r64 / 7.0_r64) ** 0.15_r64 * meteo%uw(i) !'wind speed at 10 m (m/s)
                kgnh3 = 175.287_r64 * uw10 + 262.9305_r64 !'nh3 mass transfer velocity (m/d)
                henh3 = 0.0000136785_r64 * 1.052_r64 ** (te(i, 1) - 20.0_r64) !'henry's law constant for nh3 gas (atm m^3/mole)
                vnh3 = k1nh3 * henh3 / (henh3 + 0.00008206_r64 * (te(i, 1) + 273.15_r64) &
                    * (k1nh3 / kgnh3)) !'nh3 gas transfer coefficient (m/d)
                naus = 0.000000002_r64 / henh3 * 14000.0_r64 !'sat'n conc of nh3 assuming partial press of nh3 in atm = 2e-9 atm
                !values range from 1-10e-9 in rural and 10-100e-9 in heavily polluted areas
                nh3gas = vnh3 * hydrau%reach(i)%ast * &
                    (naus - (1.0_r64 - fi) * c(i, IDX_AMMONIA_N, 1)) !'loss or gain of nh3 via gas transfer (mgn/day)

                !determine solar radiation
                iat(i) = lightheat%par * solar%jsnt(i)

                call oxygen_inhibition_and_enhancement(rates, &
                    c(i, IDX_DISSOLVED_OXYGEN, 1), fcarb, fnitr, fdenitr, frespp, frespb) !gp 27-oct-04

                !light extinction

                !gp 13-feb-06 include macrophyte extinction if bottom plants are simulated as macrophytes
                if ((hydrau%reach(i)%nupwcfrac < 1.0_r64) .or. (hydrau%reach(i)%pupwcfrac < 1.0_r64)) then
                    !include macrophyte biomass extinction if macrophytes are simulated
                    ke = lightextinction(lightheat, &
                        c(i, IDX_PARTICULATE_ORG_MAT, 1), c(i, IDX_INORG_SUSP_SOLIDS, 1), c(i, IDX_PHYTOPLANKTON, 1), &
                        c(i, IDX_BENTHIC_ALGAE, 1) / hydrau%reach(i)%depth)
                else
                    !do not include macrophyte biomass extinction of periphyton are simulated
                    ke = lightextinction(lightheat, c(i, IDX_PARTICULATE_ORG_MAT, 1), &
                        c(i, IDX_INORG_SUSP_SOLIDS, 1), c(i, IDX_PHYTOPLANKTON, 1), 0.0_r64)
                end if

                !
                ! --- (11) phytoplankton (mga/day) ---
                !

                !nutrient limitation

                !gp 09-dec-09 include fi
                !if (hydrau%reach(i)%ksn + c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1) <= 0) then
                ! phint = 0
                !else
                ! phint = (c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) / (hydrau%reach(i)%ksn + c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1))
                !end if
                if (hydrau%reach(i)%ksn + fi * c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1) <= 0) then
                    phint = 0
                else
                    phint = (fi * c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) / &
                        (hydrau%reach(i)%ksn + fi * c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1))
                end if

                if ((hydrau%reach(i)%ksp + c(i, IDX_SOLUBLE_REACTIVE_P, 1)) <= 0) then
                    phip = 0
                else
                    phip = c(i, IDX_SOLUBLE_REACTIVE_P, 1) / (hydrau%reach(i)%ksp + c(i, IDX_SOLUBLE_REACTIVE_P, 1))
                end if
                if (phip < phint) phint = phip
                if (rates%hco3use == "No") then
                    phic = alp0 * c(i, IDX_TOTAL_INORGANIC_C, 1) / (rates%ksc + alp0 * c(i, IDX_TOTAL_INORGANIC_C, 1))
                else
                    phic = (alp0 + alp1) * c(i, IDX_TOTAL_INORGANIC_C, 1) / &
                        (rates%ksc + (alp0 + alp1) * c(i, IDX_TOTAL_INORGANIC_C, 1))
                end if
                if (phic < phint) phint = phic

                !light limitation

                select case (rates%ilight)
                  case (1) !half-saturation
                    if (hydrau%reach(i)%isat == 0) then
                        phil = 1.0
                    else
                        phil = 1.0 / (ke * hydrau%reach(i)%depth) * log((hydrau%reach(i)%isat + iat(i)) &
                            / (hydrau%reach(i)%isat + iat(i) * exp(-ke * hydrau%reach(i)%depth)))
                    end if
                  case (2) !smith
                    if (hydrau%reach(i)%isat == 0) then
                        phil = 1.0
                    else
                        num = iat(i) / hydrau%reach(i)%isat + sqrt(1 + (iat(i) / hydrau%reach(i)%isat) ** 2)
                        den = iat(i) * exp(-ke * hydrau%reach(i)%depth) / hydrau%reach(i)%isat &
                            + sqrt(1 + (iat(i) * exp(-ke * hydrau%reach(i)%depth) / hydrau%reach(i)%isat) ** 2)
                        phil = 1.0 / (ke * hydrau%reach(i)%depth) * log(num / den)
                    end if
                  case (3) !steele
                    if (hydrau%reach(i)%isat == 0) then
                        phil = 1.0
                    else
                        alpha0 = iat(i) / hydrau%reach(i)%isat
                        alpha1 = iat(i) / hydrau%reach(i)%isat * exp(-ke * hydrau%reach(i)%depth)
                        phil = exp(1.0_r64) * (exp(-alpha1) - exp(-alpha0)) / (ke * hydrau%reach(i)%depth)
                    end if
                end select

                phytophoto = phil * phint * kgat(i) * hydrau%reach(i)%vol * c(i, IDX_PHYTOPLANKTON, 1) !mga/day
                phytoresp = frespp * kreat(i) * hydrau%reach(i)%vol * c(i, IDX_PHYTOPLANKTON, 1) !mga/day
                phytodeath = kdeat(i) * hydrau%reach(i)%vol * c(i, IDX_PHYTOPLANKTON, 1) !mga/day
                phytosettl = hydrau%reach(i)%va * hydrau%reach(i)%asd * c(i, IDX_PHYTOPLANKTON, 1) !mga/day
                dc(i, IDX_PHYTOPLANKTON, 1) = dc(i, IDX_PHYTOPLANKTON, 1) + phytophoto &
                    - phytoresp - phytodeath - phytosettl !mga/day

                !
                ! --- (16) bottom algae (gd/day) ---
                !

                !luxury uptake calcs

                if (c(i, IDX_BENTHIC_ALGAE, 1) > 0) then
                    ninb = inb(i) / c(i, IDX_BENTHIC_ALGAE, 1)
                    nipb = ipb(i) / c(i, IDX_BENTHIC_ALGAE, 1)
                else
                    ninb = rates%ana / rates%ada
                    nipb = rates%apa / rates%ada
                end if
                if (ninb > 0) then
                    fbnb = hydrau%reach(i)%ninbmin / ninb
                else
                    fbnb = 1.0_r64
                end if
                if (fbnb < 0) fbnb = 0.0_r64
                if (fbnb > 1.0_r64) fbnb = 1.0_r64
                if (nipb > 0.0_r64) then
                    fbpb = hydrau%reach(i)%nipbmin / nipb
                else
                    fbpb = 1.0_r64
                end if
                if (fbpb < 0) fbpb = 0.0_r64
                if (fbpb > 1.0_r64) fbpb = 1.0_r64

                phint = 1.0_r64 - fbnb
                phinsave(i) = phint !'gp 20-oct-04
                phip = 1.0_r64 - fbpb
                phipsave(i) = phip !'gp 20-oct-04
                if (phip < phint) phint = phip

                if (rates%hco3usef == "No") then
                    phic = alp0 * c(i, IDX_TOTAL_INORGANIC_C, 1) / &
                        (rates%kscf + alp0 * c(i, IDX_TOTAL_INORGANIC_C, 1)) !gp !co2 is the limiting substrate
                else
                    phic = (alp0 + alp1) * c(i, IDX_TOTAL_INORGANIC_C, 1) / &
                        (rates%kscf + (alp0 + alp1) * c(i, IDX_TOTAL_INORGANIC_C, 1)) !gp !hco3- is the limiting substrate
                end if
                phicsave(i) = phic !'gp 20-oct-04
                if (phic < phint) phint = phic

                !light limitation

                botlight = iat(i) * exp(-ke * hydrau%reach(i)%depth)
                if ((hydrau%reach(i)%nupwcfrac < 1.0_r64) .or. (hydrau%reach(i)%pupwcfrac < 1.0_r64)) then
                    !use water column average light for macrophytes
                    select case (rates%ilightf)
                      case (1) !half-saturation
                        if (hydrau%reach(i)%isatf == 0) then
                            phil = 1.0
                        else
                            phil = 1.0 / (ke * hydrau%reach(i)%depth) * log((hydrau%reach(i)%isatf + iat(i)) &
                                / (hydrau%reach(i)%isatf + iat(i) * exp(-ke * hydrau%reach(i)%depth)))
                        end if
                      case (2) !smith
                        if (hydrau%reach(i)%isatf == 0) then
                            phil = 1.0
                        else
                            num = iat(i) / hydrau%reach(i)%isatf + sqrt(1 + (iat(i) / hydrau%reach(i)%isatf) ** 2)
                            den = iat(i) * exp(-ke * hydrau%reach(i)%depth) / hydrau%reach(i)%isatf &
                                + sqrt(1 + (iat(i) * exp(-ke * hydrau%reach(i)%depth) / hydrau%reach(i)%isatf) ** 2)
                            phil = 1.0 / (ke * hydrau%reach(i)%depth) * log(num / den)
                        end if
                      case (3) !steele
                        if (hydrau%reach(i)%isatf == 0) then
                            phil = 1.0
                        else
                            alpha0 = iat(i) / hydrau%reach(i)%isatf
                            alpha1 = iat(i) / hydrau%reach(i)%isatf * exp(-ke * hydrau%reach(i)%depth)
                            phil = exp(1.0_r64) * (exp(-alpha1) - exp(-alpha0)) / (ke * hydrau%reach(i)%depth)
                        end if
                    end select
                else
                    !use bottom light for periphyton
                    select case (rates%ilightf)
                      case (1) !half-saturation
                        if (hydrau%reach(i)%isatf + botlight <= 0) then
                            phil = 1.0_r64
                        else
                            phil = botlight / (hydrau%reach(i)%isatf + botlight)
                        end if
                      case (2) !smith
                        if (botlight ** 2 + hydrau%reach(i)%isatf ** 2 == 0) then
                            phil = 1.0_r64
                        else
                            phil = botlight / sqrt(botlight ** 2 + hydrau%reach(i)%isatf ** 2)
                        end if
                      case (3) !steele
                        if (hydrau%reach(i)%isatf <= 0) then
                            phil = 1.0_r64
                        else
                            phil = botlight / hydrau%reach(i)%isatf * exp(1 - botlight / hydrau%reach(i)%isatf)
                        end if
                    end select
                end if

                philsave(i) = phil

                if (rates%typef == "Zero-order") then
                    botalgphoto = phil * phint * kgaft(i) * hydrau%reach(i)%asb !gd/day
                else
                    botalgphoto = phil * phint * kgaft(i) * hydrau%reach(i)%asb * &
                        c(i, IDX_BENTHIC_ALGAE, 1) * (1 - c(i, IDX_BENTHIC_ALGAE, 1) / hydrau%reach(i)%abmax) !gd/day
                end if

                phitotalsave(i) = phil * phint * phitsave(i)

                !basal resp
                botalgresp = krea1ft(i) * hydrau%reach(i)%asb * c(i, IDX_BENTHIC_ALGAE, 1) !gd/day
                !add phot-resp
                if (rates%typef == "Zero-order") then
                    botalgresp = botalgresp + krea2ft(i) * phil * phint * kgaft(i) * hydrau%reach(i)%asb
                else !'else first-order growth used for photo resp
                    botalgresp = botalgresp + krea2ft(i) * phil * phint * kgaft(i) * hydrau%reach(i)%asb &
                        * c(i, IDX_BENTHIC_ALGAE, 1) * (1 - c(i, IDX_BENTHIC_ALGAE, 1) / hydrau%reach(i)%abmax)
                end if
                !'adjust for oxygen attenuation
                botalgresp = frespb * botalgresp !gd/day

                botalgexc = kexaft(i) * hydrau%reach(i)%asb * c(i, IDX_BENTHIC_ALGAE, 1) !gd/day
                !note that botalgexc does not contribute to dc(i,IDX_BENTHIC_ALGAE,1)
                botalgdeath = kdeaft(i) * hydrau%reach(i)%asb * c(i, IDX_BENTHIC_ALGAE, 1) !gd/day
                dc(i, IDX_BENTHIC_ALGAE, 1) = botalgphoto - botalgresp - botalgdeath !gd/day

                !gp 25-jun-09
                savebotalgphoto(i) = botalgphoto / hydrau%reach(i)%asb
                savebotalgresp(i) = botalgresp / hydrau%reach(i)%asb
                savebotalgdeath(i) = botalgdeath / hydrau%reach(i)%asb
                savebotalgnetgrowth(i) = dc(i, IDX_BENTHIC_ALGAE, 1) / hydrau%reach(i)%asb

                !periphyton uptake of n and p

                !gp 03-dec-09
                !if ((hydrau%reach(i)%ksnf + c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) > 0) then
                ! nlim = hydrau%reach(i)%kqn / (hydrau%reach(i)%kqn + ninb - hydrau%reach(i)%ninbmin)
                ! if (nlim < 0) then
                ! nlim = 0
                ! else if (nlim > 1.0) then
                ! nlim = 1.0_r64
                ! end if
                ! botalguptaken = hydrau%reach(i)%asb * nlim * c(i, IDX_BENTHIC_ALGAE, 1) * hydrau%reach(i)%ninbupmax &
                ! * (c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) / (hydrau%reach(i)%ksnf + c(i, IDX_AMMONIA_N, 1) &
                ! + c(i, IDX_NOX_N, 1))
                !else
                ! botalguptaken = 0
                !end if
                if ((hydrau%reach(i)%ksnf + fi * c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) > 0) then
                    nlim = hydrau%reach(i)%kqn / (hydrau%reach(i)%kqn + ninb - hydrau%reach(i)%ninbmin)
                    if (nlim < 0) then
                        nlim = 0
                    else if (nlim > 1.0) then
                        nlim = 1.0_r64
                    end if
                    botalguptaken = hydrau%reach(i)%asb * nlim * c(i, IDX_BENTHIC_ALGAE, 1) * &
                        hydrau%reach(i)%ninbupmax * (fi * c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) / &
                        (hydrau%reach(i)%ksnf + fi * c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1))
                else
                    botalguptaken = 0
                end if

                if ((hydrau%reach(i)%kspf + c(i, IDX_SOLUBLE_REACTIVE_P, 1)) > 0) then
                    plim = hydrau%reach(i)%kqp / (hydrau%reach(i)%kqp + nipb - hydrau%reach(i)%nipbmin)
                    if (plim < 0) then
                        plim = 0
                    else if (plim > 1) then
                        plim = 1.0_r64
                    end if
                    botalguptakep = hydrau%reach(i)%asb * plim * c(i, IDX_BENTHIC_ALGAE, 1) * &
                        hydrau%reach(i)%nipbupmax * c(i, IDX_SOLUBLE_REACTIVE_P, 1) / &
                        (hydrau%reach(i)%kspf + c(i, IDX_SOLUBLE_REACTIVE_P, 1))
                else
                    botalguptakep = 0
                end if

                !gp 03-dec-09
                !if (botalguptaken * sys%dt > (c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) * hydrau%reach(i)%vol) then
                ! botalguptaken = (c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) * hydrau%reach(i)%vol / sys%dt !'mgn/day
                !end if
                !if (botalguptakep * sys%dt > c(i, IDX_SOLUBLE_REACTIVE_P, 1) * hydrau%reach(i)%vol) then
                ! botalguptakep = c(i, IDX_SOLUBLE_REACTIVE_P, 1) * hydrau%reach(i)%vol / sys%dt !'mgp/day
                !end if
                if (botalguptaken * sys%dt > (fi * c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) * hydrau%reach(i)%vol) then
                    botalguptaken = (fi * c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) * hydrau%reach(i)%vol / sys%dt !'mgn/day
                end if
                if (botalguptakep * sys%dt > c(i, IDX_SOLUBLE_REACTIVE_P, 1) * hydrau%reach(i)%vol) then
                    botalguptakep = c(i, IDX_SOLUBLE_REACTIVE_P, 1) * hydrau%reach(i)%vol / sys%dt !'mgp/day
                end if

                !change in intracellular n and p in periphyton
                dinb(i) = botalguptaken - ninb * botalgdeath - ninb * botalgexc !mgn/day
                dipb(i) = botalguptakep - nipb * botalgdeath - nipb * botalgexc !mgp/day

                !ammonium preference

                !gp 03-dec-09
                !prefam = 0
                !if (c(i, IDX_AMMONIA_N, 1) * c(i, IDX_NOX_N, 1) > 0) then
                ! prefam = c(i, IDX_AMMONIA_N, 1) * c(i, IDX_NOX_N, 1) / (hydrau%reach(i)%khnx + c(i, IDX_AMMONIA_N, 1)) / (hydrau%reach(i)%khnx + c(i, IDX_NOX_N, 1)) &
                ! + c(i, IDX_AMMONIA_N, 1) * hydrau%reach(i)%khnx / (c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) / (hydrau%reach(i)%khnx + c(i, IDX_NOX_N, 1))
                !end if
                !prefamf = 0
                !if (c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1) > 0) then
                ! prefamf = c(i, IDX_AMMONIA_N, 1) * c(i, IDX_NOX_N, 1) / (hydrau%reach(i)%khnxf + c(i, IDX_AMMONIA_N, 1)) / (hydrau%reach(i)%khnxf + c(i, IDX_NOX_N, 1)) &
                ! + c(i, IDX_AMMONIA_N, 1) * hydrau%reach(i)%khnxf / (c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) / (hydrau%reach(i)%khnxf + c(i, IDX_NOX_N, 1))
                !end if
                prefam = 0
                if (fi * c(i, IDX_AMMONIA_N, 1) * c(i, IDX_NOX_N, 1) > 0) then
                    prefam = fi * c(i, IDX_AMMONIA_N, 1) * c(i, IDX_NOX_N, 1) / &
                        (hydrau%reach(i)%khnx + fi * c(i, IDX_AMMONIA_N, 1)) / &
                        (hydrau%reach(i)%khnx + c(i, IDX_NOX_N, 1)) &
                        + fi * c(i, IDX_AMMONIA_N, 1) * hydrau%reach(i)%khnx / &
                        (fi * c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) / &
                        (hydrau%reach(i)%khnx + c(i, IDX_NOX_N, 1))
                end if
                prefamf = 0
                if (fi * c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1) > 0) then
                    prefamf = fi * c(i, IDX_AMMONIA_N, 1) * c(i, IDX_NOX_N, 1) / &
                        (hydrau%reach(i)%khnxf + fi * c(i, IDX_AMMONIA_N, 1)) / &
                        (hydrau%reach(i)%khnxf + c(i, IDX_NOX_N, 1)) &
                        + fi * c(i, IDX_AMMONIA_N, 1) * hydrau%reach(i)%khnxf / &
                        (fi * c(i, IDX_AMMONIA_N, 1) + c(i, IDX_NOX_N, 1)) / &
                        (hydrau%reach(i)%khnxf + c(i, IDX_NOX_N, 1))
                end if

                !
                ! --- sediment fluxes of o2, c, n, and p ---
                !

                if (sys%calcsedflux == "Yes" .or. sys%calcsedflux == "Option 1" .or. &
                    sys%calcsedflux == "Option 2") then !gp 11-jan-06

                    Jcin = Rates%roc * Rates%aca * hydrau%reach(i)%va * c(i, IDX_PHYTOPLANKTON, 1) &
                        + hydrau%reach(i)%vdt * c(i, IDX_PARTICULATE_ORG_MAT, 1) / Rates%adc * Rates%roc !gO2/m^2/d
                    Jnin = (Rates%ana * hydrau%reach(i)%va * c(i, IDX_PHYTOPLANKTON, 1) + &
                        hydrau%reach(i)%von * c(i, IDX_ORGANIC_N, 1)) / 1000.0_r64 !gN/m^2/d
                    Jpin = (Rates%apa * hydrau%reach(i)%va * c(i, IDX_PHYTOPLANKTON, 1) + &
                        hydrau%reach(i)%vop * c(i, IDX_ORGANIC_P, 1) + hydrau%reach(i)%vip * c(i, IDX_SOLUBLE_REACTIVE_P, 1)) / &
                        1000.0_r64 !gP/m^2/d
                    ow = c(i, IDX_DISSOLVED_OXYGEN, 1) !mgO2/L = gO2/m^3

                    !gp 09-Dec-09 include Fi
                    !NH3w = c(i, IDX_AMMONIA_N, 1) / 1000.0_r64 !gN/m^3
                    NH3w = Fi * c(i, IDX_AMMONIA_N, 1) / 1000.0_r64 !gN/m^3

                    NO3w = c(i, IDX_NOX_N, 1) / 1000.0_r64 !gN/m^3
                    PO4w = c(i, IDX_SOLUBLE_REACTIVE_P, 1) / 1000.0_r64 !gP/m^3
                    Tw = Te(i, 1) !deg C

                    !note that inputs and outputs of O2, C(O2), N, and P fluxes to/from SedCalcNumNew are g/m^2/day
                    CALL SedCalcNumNew(Jcin, Jnin, Jpin, ow, hydrau%reach(i)%depth, Tw, SOD, Jamm, Jnitr, Jmeth, Jmethg, &
                        Jphos, NH3w, NO3w, PO4w, c(i, IDX_CBOD_FAST, 1), CSOD, sys%calcSedFlux)

                    JNH4 = Jamm * 1000.0_r64 !mgN/m^2/d
                    JNO3 = Jnitr * 1000.0_r64 !mgN/m^2/d
                    JCH4 = Jmeth !gO2/m^2/d
                    JPO4 = Jphos * 1000.0_r64 * Rates%kspi / (Rates%kspi + c(i, IDX_DISSOLVED_OXYGEN, 1)) !mgP/m^2/d
                END IF

                SODpr(i) = hydrau%reach(i)%SODspec + SOD !gO2/m^2/d
                JNH4pr(i) = hydrau%reach(i)%JNH4spec + JNH4 !mgN/m^2/d
                JNO3pr(i) = JNO3 !mgN/m^2/d
                JCH4pr(i) = hydrau%reach(i)%JCH4spec + JCH4 !gO2/m^2/d
                JSRPpr(i) = hydrau%reach(i)%JSRPspec + JPO4 !mgP/m^2/d
                CSODpr(i) = CSOD !gO2/m^2/d

                !gp 01-Nov-04 output diel diagenesis fluxes of constituents (averaged over total surface area)
                DiagFluxDO(i) = -SODpr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast !'dissolved oxygen gO2/m^2/d
                DiagFluxCBOD(i) = JCH4pr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast !'fast CBOD gO2/m^2/d
                DiagFluxNH4(i) = JNH4pr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast !'ammonia mgN/m^2/d
                DiagFluxNO3(i) = JNO3pr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast !'nitrate mgN/m^2/d
                DiagFluxSRP(i) = JSRPpr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast !'SRP mgP/m^2/d
                DiagFluxIC(i) = CSODpr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast / Rates%roc !'inorganic C gC/m^2/d

                !
                ! --- (12) detritus (gD/day) ---
                !

                DetrDiss = kdtT(i) * hydrau%reach(i)%vol * c(i, IDX_PARTICULATE_ORG_MAT, 1)
                DetrSettl = hydrau%reach(i)%vdt * hydrau%reach(i)%Asd * c(i, IDX_PARTICULATE_ORG_MAT, 1)
                dc(i, IDX_PARTICULATE_ORG_MAT, 1) = dc(i, IDX_PARTICULATE_ORG_MAT, 1) - DetrDiss
                dc(i, IDX_PARTICULATE_ORG_MAT, 1) = dc(i, IDX_PARTICULATE_ORG_MAT, 1) - DetrSettl
                dc(i, IDX_PARTICULATE_ORG_MAT, 1) = dc(i, IDX_PARTICULATE_ORG_MAT, 1) + Rates%ada * PhytoDeath
                dc(i, IDX_PARTICULATE_ORG_MAT, 1) = dc(i, IDX_PARTICULATE_ORG_MAT, 1) + BotAlgDeath

                !
                ! --- (6) Organic Nitrogen (mgN/day) ---
                !

                OrgNHydr = khnT(i) * hydrau%reach(i)%vol * c(i, IDX_ORGANIC_N, 1)
                OrgNSettl = hydrau%reach(i)%von * hydrau%reach(i)%Asd * c(i, IDX_ORGANIC_N, 1)
                dc(i, IDX_ORGANIC_N, 1) = dc(i, IDX_ORGANIC_N, 1) + NINb * BotAlgDeath + Rates%ana * PhytoDeath
                dc(i, IDX_ORGANIC_N, 1) = dc(i, IDX_ORGANIC_N, 1) - OrgNHydr - OrgNSettl

                !
                ! --- (7) Total Ammonia (NH4+ + NH3) Nitrogen (mgN/day) ---
                !

                !gp 03-Dec-09
                !NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%vol * c(i, IDX_AMMONIA_N, 1)
                NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%vol * Fi * c(i, IDX_AMMONIA_N, 1)

                dc(i, IDX_AMMONIA_N, 1) = dc(i, IDX_AMMONIA_N, 1) + OrgNHydr
                dc(i, IDX_AMMONIA_N, 1) = dc(i, IDX_AMMONIA_N, 1) - NH4Nitrif
                dc(i, IDX_AMMONIA_N, 1) = dc(i, IDX_AMMONIA_N, 1) + Rates%ana * PhytoResp
                dc(i, IDX_AMMONIA_N, 1) = dc(i, IDX_AMMONIA_N, 1) - prefam * Rates%ana * PhytoPhoto
                dc(i, IDX_AMMONIA_N, 1) = dc(i, IDX_AMMONIA_N, 1) - prefamF * BotAlgUptakeN * hydrau%reach(i)%NUpWCfrac
                dc(i, IDX_AMMONIA_N, 1) = dc(i, IDX_AMMONIA_N, 1) + JNH4pr(i) * hydrau%reach(i)%Asd
                dc(i, IDX_AMMONIA_N, 1) = dc(i, IDX_AMMONIA_N, 1) + NINb * BotAlgExc

                !gp 03-Dec-09
                dc(i, IDX_AMMONIA_N, 1) = dc(i, IDX_AMMONIA_N, 1) + NH3gas !air/water exchange of NH3 gas

                !
                ! --- (8) Nitrate+Nitrite Nitrogen (mgN/day) ---
                !

                Denitr = c(i, IDX_CBOD_FAST, 1) / (0.1_r64 + c(i, IDX_CBOD_FAST, 1)) * fdenitr * &
                    kiT(i) * hydrau%reach(i)%vol * c(i, IDX_NOX_N, 1) !wc denitr
                dc(i, IDX_NOX_N, 1) = dc(i, IDX_NOX_N, 1) + NH4Nitrif
                dc(i, IDX_NOX_N, 1) = dc(i, IDX_NOX_N, 1) - Denitr
                dc(i, IDX_NOX_N, 1) = dc(i, IDX_NOX_N, 1) - (1 - prefam) * Rates%ana * PhytoPhoto
                dc(i, IDX_NOX_N, 1) = dc(i, IDX_NOX_N, 1) - (1 - prefamF) * BotAlgUptakeN * hydrau%reach(i)%NUpWCfrac
                dc(i, IDX_NOX_N, 1) = dc(i, IDX_NOX_N, 1) + JNO3pr(i) * hydrau%reach(i)%Asd
                dc(i, IDX_NOX_N, 1) = dc(i, IDX_NOX_N, 1) - vdiT(i) * hydrau%reach(i)%Asd * c(i, IDX_NOX_N, 1) !sed denitr

                !
                ! --- (9) Organic Phosphorus (mgP/day) ---
                !

                OrgPHydr = khpT(i) * hydrau%reach(i)%vol * c(i, IDX_ORGANIC_P, 1)
                OrgPSettl = hydrau%reach(i)%vop * hydrau%reach(i)%Asd * c(i, IDX_ORGANIC_P, 1)
                dc(i, IDX_ORGANIC_P, 1) = dc(i, IDX_ORGANIC_P, 1) + NIPb * BotAlgDeath + Rates%apa * PhytoDeath
                dc(i, IDX_ORGANIC_P, 1) = dc(i, IDX_ORGANIC_P, 1) - OrgPHydr - OrgPSettl

                !
                ! --- (10) Inorganic Soluble Reactive Phosphorus (mgP/day) ---
                !

                dc(i, IDX_SOLUBLE_REACTIVE_P, 1) = dc(i, IDX_SOLUBLE_REACTIVE_P, 1) + OrgPHydr
                dc(i, IDX_SOLUBLE_REACTIVE_P, 1) = dc(i, IDX_SOLUBLE_REACTIVE_P, 1) + Rates%apa * PhytoResp
                dc(i, IDX_SOLUBLE_REACTIVE_P, 1) = dc(i, IDX_SOLUBLE_REACTIVE_P, 1) - Rates%apa * PhytoPhoto
                dc(i, IDX_SOLUBLE_REACTIVE_P, 1) = dc(i, IDX_SOLUBLE_REACTIVE_P, 1) - BotAlgUptakeP * hydrau%reach(i)%PUpWCfrac
                dc(i, IDX_SOLUBLE_REACTIVE_P, 1) = dc(i, IDX_SOLUBLE_REACTIVE_P, 1) + JSRPpr(i) * hydrau%reach(i)%Asd
                InorgPSettl = hydrau%reach(i)%vip * hydrau%reach(i)%Asd * c(i, IDX_SOLUBLE_REACTIVE_P, 1)
                dc(i, IDX_SOLUBLE_REACTIVE_P, 1) = dc(i, IDX_SOLUBLE_REACTIVE_P, 1) - InorgPSettl
                dc(i, IDX_SOLUBLE_REACTIVE_P, 1) = dc(i, IDX_SOLUBLE_REACTIVE_P, 1) + NIPb * BotAlgExc

                !
                ! --- (4) CBOD Slow (gO2/day) ---
                !

                CBODsHydr = khcT(i) * hydrau%reach(i)%vol * c(i, IDX_CBOD_SLOW, 1)
                dc(i, IDX_CBOD_SLOW, 1) = dc(i, IDX_CBOD_SLOW, 1) + Rates%roc * (1.0_r64 / Rates%adc) * DetrDiss
                dc(i, IDX_CBOD_SLOW, 1) = dc(i, IDX_CBOD_SLOW, 1) - CBODsHydr
                CBODsOxid = fcarb * kdcsT(i) * hydrau%reach(i)%vol * c(i,IDX_CBOD_SLOW, 1)
                dc(i, IDX_CBOD_SLOW, 1) = dc(i, IDX_CBOD_SLOW, 1) - CBODsOxid

                !
                ! --- (5) CBOD Fast (gO2/day) ---
                !

                CBODfOxid = fcarb * kdcT(i) * hydrau%reach(i)%vol * c(i, IDX_CBOD_FAST, 1)
                dc(i, IDX_CBOD_FAST, 1) = dc(i, IDX_CBOD_FAST, 1) + CBODsHydr
                dc(i, IDX_CBOD_FAST, 1) = dc(i, IDX_CBOD_FAST, 1) - CBODfOxid
                dc(i, IDX_CBOD_FAST, 1) = dc(i, IDX_CBOD_FAST, 1) + JCH4pr(i) * hydrau%reach(i)%Asd
                dc(i, IDX_CBOD_FAST, 1) = dc(i, IDX_CBOD_FAST, 1) - Rates%rondn * Denitr !gp end new block

                !
                ! --- (2) Inorganic Suspended Solids (gD/day) ---
                !

                InorgSettl = hydrau%reach(i)%vss * hydrau%reach(i)%Asd * c(i, IDX_INORG_SUSP_SOLIDS, 1)
                dc(i, IDX_INORG_SUSP_SOLIDS, 1) = dc(i, IDX_INORG_SUSP_SOLIDS, 1) - InorgSettl !gp end new block

                !
                ! --- (14) Generic constituent or COD (user defined units of mass/time, gO2/day if used as COD) ---
                !

                dc(i, IDX_GENERIC_CONST, 1) = dc(i, IDX_GENERIC_CONST, 1) - kgenT(i) * &
                    hydrau%reach(i)%vol * c(i, IDX_GENERIC_CONST, 1)
                dc(i, IDX_GENERIC_CONST, 1) = dc(i, IDX_GENERIC_CONST, 1) - hydrau%reach(i)%vgen * &
                    hydrau%reach(i)%Asd * c(i, IDX_GENERIC_CONST, 1)
                IF (Rates%useGenericAsCOD == "Yes") THEN
                    CODoxid = kgenT(i) * hydrau%reach(i)%vol * c(i, IDX_GENERIC_CONST, 1)
                ELSE
                    CODoxid = 0
                END IF

                !
                ! --- (3) Dissolved Oxygen (gO2/day) ---
                !

                OxReaer = kaT(i) * hydrau%reach(i)%vol * &
                    (oxygen_saturation(Te(i, 1), hydrau%reach(i)%elev) - c(i, IDX_DISSOLVED_OXYGEN, 1))
                dc(i, IDX_DISSOLVED_OXYGEN, 1) = dc(i, IDX_DISSOLVED_OXYGEN, 1) + OxReaer
                dc(i, IDX_DISSOLVED_OXYGEN, 1) = dc(i, IDX_DISSOLVED_OXYGEN, 1) - CBODsOxid - CBODfOxid
                dc(i, IDX_DISSOLVED_OXYGEN, 1) = dc(i, IDX_DISSOLVED_OXYGEN, 1) - Rates%ron * NH4Nitrif
                dc(i, IDX_DISSOLVED_OXYGEN, 1) = dc(i, IDX_DISSOLVED_OXYGEN, 1) - Rates%roa * PhytoResp
                dc(i, IDX_DISSOLVED_OXYGEN, 1) = dc(i, IDX_DISSOLVED_OXYGEN, 1) + Rates%roa * PhytoPhoto * prefam
                dc(i, IDX_DISSOLVED_OXYGEN, 1) = dc(i, IDX_DISSOLVED_OXYGEN, 1) + &
                    Rates%roa * PhytoPhoto * (1.0_r64 - prefam) * 138.0_r64 / 107.0_r64
                dc(i, IDX_DISSOLVED_OXYGEN, 1) = dc(i, IDX_DISSOLVED_OXYGEN, 1) - Rates%roc / Rates%adc * BotAlgResp
                dc(i, IDX_DISSOLVED_OXYGEN, 1) = dc(i, IDX_DISSOLVED_OXYGEN, 1) + Rates%roc / Rates%adc * BotAlgPhoto * prefamF
                dc(i, IDX_DISSOLVED_OXYGEN, 1) = dc(i, IDX_DISSOLVED_OXYGEN, 1) + &
                    Rates%roc / Rates%adc * BotAlgPhoto * (1 - prefamF) * 138.0_r64 / 107.0_r64
                dc(i, IDX_DISSOLVED_OXYGEN, 1) = dc(i, IDX_DISSOLVED_OXYGEN, 1) - SODpr(i) * hydrau%reach(i)%Asd
                dc(i, IDX_DISSOLVED_OXYGEN, 1) = dc(i, IDX_DISSOLVED_OXYGEN, 1) - CODoxid

                !'gp 05-Jul-05 save DO fluxes (gO2/m^2/d)
                saveDOfluxReaer(i) = OxReaer / hydrau%reach(i)%Ast
                saveDOfluxCBODfast(i) = -CBODfOxid / hydrau%reach(i)%Ast
                saveDOfluxCBODslow(i) = -CBODsOxid / hydrau%reach(i)%Ast
                saveDOfluxNitrif(i) = -Rates%ron * NH4Nitrif / hydrau%reach(i)%Ast
                saveDOfluxPhytoResp(i) = -Rates%roa * PhytoResp / hydrau%reach(i)%Ast
                saveDOfluxPhytoPhoto(i) = (Rates%roa * PhytoPhoto * prefam &
                    + Rates%roa * PhytoPhoto * (1 - prefam) * 138.0_r64 / 107.0_r64) &
                    / hydrau%reach(i)%Ast
                saveDOfluxBotalgResp(i) = -Rates%roc / Rates%adc * BotAlgResp / hydrau%reach(i)%Ast
                saveDOfluxBotalgPhoto(i) = (Rates%roc / Rates%adc * BotAlgPhoto * prefamF &
                    + Rates%roc / Rates%adc * BotAlgPhoto * (1 - prefamF) * 138.0_r64 / 107.0_r64) &
                    / hydrau%reach(i)%Ast
                saveDOfluxSOD(i) = -SODpr(i) * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast
                saveDOfluxCOD(i) = -CODoxid / hydrau%reach(i)%Ast

                dc(i, IDX_TOTAL_INORGANIC_C, 1) = dc(i, IDX_TOTAL_INORGANIC_C, 1) - Rates%rcca * PhytoPhoto
                dc(i, IDX_TOTAL_INORGANIC_C, 1) = dc(i, IDX_TOTAL_INORGANIC_C, 1) + Rates%rcca * PhytoResp
                dc(i, IDX_TOTAL_INORGANIC_C, 1) = dc(i, IDX_TOTAL_INORGANIC_C, 1) + Rates%rcco * CBODfOxid
                dc(i, IDX_TOTAL_INORGANIC_C, 1) = dc(i, IDX_TOTAL_INORGANIC_C, 1) + Rates%rcco * CBODsOxid
                dc(i, IDX_TOTAL_INORGANIC_C, 1) = dc(i, IDX_TOTAL_INORGANIC_C, 1) - Rates%rccd * BotAlgPhoto
                dc(i, IDX_TOTAL_INORGANIC_C, 1) = dc(i, IDX_TOTAL_INORGANIC_C, 1) + Rates%rccd * BotAlgResp
                dc(i, IDX_TOTAL_INORGANIC_C, 1) = dc(i, IDX_TOTAL_INORGANIC_C, 1) + Rates%rcco * CSOD * hydrau%reach(i)%Asd
                dc(i, IDX_TOTAL_INORGANIC_C, 1) = dc(i, IDX_TOTAL_INORGANIC_C, 1) + kacT(i) * hydrau%reach(i)%vol * (Khs(i, 1) &
                    * Rates%pco2 - alp0 * c(i, IDX_TOTAL_INORGANIC_C, 1)) !gp end new block

                !
                ! --- (nv-2) Alkalinity (gCaCO3/day) ---
                !

                If (sys%simAlk == "Yes") Then

                    !gp 03-Dec-09
                    !dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkaa * PhytoPhoto * prefam * 50000.0_r64
                    !dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkan * PhytoPhoto * (1.0_r64 - prefam) * 50000.0_r64
                    !dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkaa * PhytoResp * 50000.0_r64
                    !dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkbn * BotAlgUptakeN * prefamF * 50000.0_r64 &
                    ! - Rates%ralkbp * BotAlgUptakeP * 50000.0_r64
                    !dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbn * BotAlgUptakeN * (1.0_r64 - prefamF) * 50000.0_r64
                    !dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkn * NH4Nitrif * 50000.0_r64
                    !dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkden * Denitr * 50000.0_r64
                    !dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbn * OrgNHydr * 50000.0_r64 + Rates%ralkbp * OrgPHydr * 50000.0_r64 !gp end new block
                    !dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkn * JNH4pr(i) * hydrau%reach(i)%Asd * 50000.0_r64 !'sed flux of ammonia
                    !dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkbp * JSRPpr(i) * hydrau%reach(i)%Asd * 50000.0_r64 !'sed flux of PO4

                    !gp 03-Dec-09
                    ! alkalinity derivative (factor of 50043.45 converts eqH+/L to mgCaCO3/L or gCaCO3/m^3)
                    ! Fi is already accounted for in nitrification, uptake
                    ! but not in the production of total ammonia (e.g. excretion, OrgN hydrolysis, JNH4pr)
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkbn * Fi * Rates%ana * &
                        PhytoPhoto * prefam * 50043.45_r64 !'phyto photo uptake of ammonia
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbn * Rates%ana * &
                        PhytoPhoto * (1 - prefam) * 50043.45_r64 !'phyto photo uptake of nitrate
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbp * Pcharge * &
                        Rates%apa * PhytoPhoto * 50043.45_r64 !'phyto photo uptake of P
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbn * Fi * Rates%ana * &
                        PhytoResp * 50043.45_r64 !'phyto resp of N
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkbp * Pcharge * &
                        Rates%apa * PhytoResp * 50043.45_r64 !'phyto resp of P
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkbn * Fi * &
                        BotAlgUptakeN * prefamF * hydrau%reach(i)%NUpWCfrac * 50043.45_r64 !'periphyton N uptake (ammonia)
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbn * BotAlgUptakeN * &
                        (1 - prefamF) * hydrau%reach(i)%NUpWCfrac * 50043.45_r64 !'periphyton N uptake (nitrate)
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbp * Pcharge * &
                        BotAlgUptakeP * hydrau%reach(i)%PUpWCfrac * 50043.45_r64 !'periphyton P uptake
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbn * Fi * &
                        BotAlgExc * NINb * 50043.45_r64 !'periphyton excretion of N
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkbp * Pcharge * &
                        BotAlgExc * NIPb * 50043.45_r64 !'periphyton excretion of P
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkbn * 2.0_r64 * &
                        NH4Nitrif * 50043.45_r64 !'nitrification ammonia loss and nitrate gain
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbn * &
                        Denitr * 50043.45_r64 !'water column denitrification
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbn * vdiT(i) * &
                        hydrau%reach(i)%Asd * c(i, IDX_NOX_N, 1) * 50043.45_r64 !'sediment denitrification
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbn * Fi * &
                        OrgNHydr * 50043.45_r64 !'organic N hydrolysis
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkbp * Pcharge * &
                        OrgPHydr * 50043.45_r64 !'organic P hydrolysis
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbn * Fi * &
                        JNH4pr(i) * hydrau%reach(i)%Asd * 50043.45_r64 !'sed flux of ammonia
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkbn * &
                        JNO3pr(i) * hydrau%reach(i)%Asd * 50043.45_r64 !'sed flux of nitrate
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) - Rates%ralkbp * Pcharge * &
                        JSRPpr(i) * hydrau%reach(i)%Asd * 50043.45_r64 !'sed flux of PO4
                    dc(i, IDX_ALKALINITY, 1) = dc(i, IDX_ALKALINITY, 1) + Rates%ralkbp * Pcharge * &
                        InorgPSettl * 50043.45_r64 !'settling flux of PO4

                end if

                !'gp 05-Jul-05 save CO2 fluxes (gC/m^2/d)
                saveCO2fluxReaer(i) = (kacT(i) * hydrau%reach(i)%vol &
                    * (Khs(i, 1) * Rates%pco2 - alp0 * c(i, IDX_TOTAL_INORGANIC_C, 1))) &
                    / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxCBODfast(i) = Rates%rcco * CBODfOxid / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxCBODslow(i) = Rates%rcco * CBODsOxid / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxPhytoResp(i) = Rates%rcca * PhytoResp / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxPhytoPhoto(i) = -Rates%rcca * PhytoPhoto / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxBotalgResp(i) = Rates%rccd * BotAlgResp / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxBotalgPhoto(i) = -Rates%rccd * BotAlgPhoto / hydrau%reach(i)%Ast / Rates%rccc
                saveCO2fluxSOD(i) = Rates%rcco * CSOD * hydrau%reach(i)%Asd / hydrau%reach(i)%Ast / Rates%rccc

                !
                ! --- Pathogen indicator bacteria (13) ---
                !

                ksol = hydrau%reach(i)%apath * Solar%Jsnt(i) / 24.0_r64 &
                    * (1.0_r64 - EXP(-ke * hydrau%reach(i)%depth)) / (ke * hydrau%reach(i)%depth)
                dc(i, IDX_PATHOGEN, 1) = dc(i, IDX_PATHOGEN, 1) - (kpathT(i) + ksol) * hydrau%reach(i)%vol * c(i, IDX_PATHOGEN, 1)
                dc(i, IDX_PATHOGEN, 1) = dc(i, IDX_PATHOGEN, 1) - hydrau%reach(i)%vpath * &
                    hydrau%reach(i)%Asd * c(i, IDX_PATHOGEN, 1)

            END DO


            !
            ! -------------------------------------------------
            ! -------------------------------------------------
            ! -------------------------------------------------
            ! --- Hyporheic pore water quality deriviatives ---
            ! -------------------------------------------------
            ! -------------------------------------------------
            ! -------------------------------------------------
            !

            SELECT CASE (sys%simHyporheicWQ)


                !
                ! -----------------------------------------------------
                ! -----------------------------------------------------
                ! --- Level 1 growth kinetics for hyporheic biofilm ---
                ! -----------------------------------------------------
                ! -----------------------------------------------------
                !

              CASE ('Level 1') !gp 15-Nov-04

                !
                ! ------------------------
                ! --- calc and save pH ---
                ! ------------------------
                !

                DO i=0, nr

                    call ph_solver(sys%imethph, ph, c(i, IDX_TOTAL_INORGANIC_C, 2), Te(i, 2), &
                        c(i, IDX_ALKALINITY, 2), c(i, IDX_CONDUCTIVITY, 2))

                    pHs(i, 2) = pH
                    CALL ChemRates(Te(i, 2), K1, K2, KW, Kh, c(i, IDX_CONDUCTIVITY, 2))
                    K1s(i, 2) = K1; K2s(i, 2) = K2; Khs(i, 2) = Kh
                END DO

                !
                ! -----------------------------------------------------
                ! --- temperature adjustment of rates and constants ---
                ! -----------------------------------------------------
                !

                DO i = 1, nr

                    !gp 08-Feb-06
                    !khcT(i) = Rates%khc * Rates%tkhc ** (Te(i, 2) - 20)
                    !kdcsT(i) = Rates%kdcs * Rates%tkdcs ** (Te(i, 2) - 20)
                    !kdcT(i) = Rates%kdc * Rates%tkdc ** (Te(i, 2) - 20)
                    !khnT(i) = Rates%khn * Rates%tkhn ** (Te(i, 2) - 20)
                    !khpT(i) = Rates%khp * Rates%tkhp ** (Te(i, 2) - 20)
                    !knT(i) = Rates%kn * Rates%tkn ** (Te(i, 2) - 20)
                    !kiT(i) = Rates%ki * Rates%tki ** (Te(i, 2) - 20)
                    !os(i, 2)=oxsat(Te(i, 2), hydrau%reach(i)%elev)
                    !kdeaT(i) = Rates%kdea * Rates%tkdea ** (Te(i, 2) - 20)
                    !kreaT(i) = Rates%krea * Rates%tkrea ** (Te(i, 2) - 20)
                    !kdtT(i) = Rates%kdt * Rates%tkdt ** (Te(i, 2) - 20)
                    !kpathT(i) = Rates%kpath * Rates%tkpath ** (Te(i, 2) - 20)
                    !kgaHT(i) = Rates%kgaH * Rates%tkgaH ** (Te(i, 2) - 20) !enhanced CBOD oxidation rate from hyporheic biofilm
                    !kgenT(i) = Rates%kgen * Rates%tkgen ** (Te(i, 2) - 20) !gp 30-Nov-04
                    khcT(i) = hydrau%reach(i)%khc * Rates%tkhc ** (Te(i, 2) - 20)
                    kdcsT(i) = hydrau%reach(i)%kdcs * Rates%tkdcs ** (Te(i, 2) - 20)
                    kdcT(i) = hydrau%reach(i)%kdc * Rates%tkdc ** (Te(i, 2) - 20)
                    khnT(i) = hydrau%reach(i)%khn * Rates%tkhn ** (Te(i, 2) - 20)
                    khpT(i) = hydrau%reach(i)%khp * Rates%tkhp ** (Te(i, 2) - 20)
                    knT(i) = hydrau%reach(i)%kn * Rates%tkn ** (Te(i, 2) - 20)
                    kiT(i) = hydrau%reach(i)%ki * Rates%tki ** (Te(i, 2) - 20)
                    os(i, 2)=oxygen_saturation(Te(i, 2), hydrau%reach(i)%elev)
                    kdeaT(i) = hydrau%reach(i)%kdea * Rates%tkdea ** (Te(i, 2) - 20)
                    kreaT(i) = hydrau%reach(i)%krea * Rates%tkrea ** (Te(i, 2) - 20)
                    kdtT(i) = hydrau%reach(i)%kdt * Rates%tkdt ** (Te(i, 2) - 20)
                    kpathT(i) = hydrau%reach(i)%kpath * Rates%tkpath ** (Te(i, 2) - 20)
                    kgaHT(i) = hydrau%reach(i)%kgaH * Rates%tkgaH ** (Te(i, 2) - 20)
                    kgenT(i) = hydrau%reach(i)%kgen * Rates%tkgen ** (Te(i, 2) - 20)

                END DO

                !
                ! ---------------------------------------------------------------------------------------------
                ! --- mass transfer between water column and hyporheic zone by bulk hyporheic exchange flow ---
                ! ---------------------------------------------------------------------------------------------
                !

                DO i = 1, nr
                    DO k = 1, nv-1
                        !gp un-comment next two lines to activate mass transfer between water column and hyporheic pore water
                        dc(i, k, 1) = dc(i, k, 1) + hydrau%reach(i)%EhyporheicCMD * (c(i, k, 2) - c(i, k, 1))
                        dc(i, k, 2) = hydrau%reach(i)%EhyporheicCMD * (c(i, k, 1) - c(i, k, 2))
                    END DO
                    !store hyporheic fluxes for output (positive flux is source to the water column from the hyporheic zone)
                    HypoFluxDO(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_DISSOLVED_OXYGEN, 2) - c(i, IDX_DISSOLVED_OXYGEN, 1)) / hydrau%reach(i)%Ast
                    ! DO gO2/m^2/d
                    HypoFluxCBOD(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_CBOD_FAST, 2) - c(i, IDX_CBOD_FAST, 1)) / hydrau%reach(i)%Ast
                    ! fast CBOD gO2/m^2/d
                    HypoFluxNH4(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_AMMONIA_N, 2) - c(i, IDX_AMMONIA_N, 1)) / hydrau%reach(i)%Ast
                    ! NH4 mgN/m^2/d
                    HypoFluxNO3(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_NOX_N, 2) - c(i, IDX_NOX_N, 1)) / hydrau%reach(i)%Ast
                    ! NO3 mgN/m^2/d
                    HypoFluxSRP(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_SOLUBLE_REACTIVE_P, 2) - c(i, IDX_SOLUBLE_REACTIVE_P, 1)) / hydrau%reach(i)%Ast
                    ! SRP mgP/m^2/d
                    HypoFluxIC(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_TOTAL_INORGANIC_C, 2) - c(i, IDX_TOTAL_INORGANIC_C, 1)) / &
                        hydrau%reach(i)%Ast / Rates%rccc !cT gC/m^2/d
                END DO

                !
                ! -----------------------------------------------------
                ! --- mass derivatives for the hyporheic pore water ---
                ! -----------------------------------------------------
                !

                DO i = 1, nr

                    !
                    ! --- pH dependent variables for later use in deriv calcs in this i loop
                    !
                    !'gp 03-Dec-09
                    !'fraction of ionized ammonia
                    hh = 10.0_r64 ** (-pHs(i, 2))
                    Kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (Te(i, 2) + 273.15_r64)))
                    Fi = hh / (hh + Kamm)
                    !'phosphate fractions of H2PO4- (FPO41), HPO4-- (FPO42), and PO4--- (FPO43)
                    KPO41 = 10.0_r64 ** (-2.15_r64)
                    KPO42 = 10.0_r64 ** (-7.2_r64)
                    KPO43 = 10.0_r64 ** (-12.35_r64)
                    DPO = 1.0_r64 / (hh ** 3.0_r64 + KPO41 * hh ** 2.0_r64 + KPO41 * KPO42 * hh + KPO41 * KPO42 * KPO43)
                    FPO41 = KPO41 * hh ** 2.0_r64 * DPO !'fraction of phosphate as H2PO4-
                    FPO42 = KPO41 * KPO42 * hh * DPO !'fraction of phosphate as HPO4--
                    FPO43 = KPO41 * KPO42 * KPO43 * DPO !'fraction of phosphate as PO4---
                    !'stoichiometric conversion factors (some are pH dependent)
                    !ralkbn = 1# / 14.0067 / 1000# / 1000# !'eqH+/L per mgN/m^3
                    !ralkbp = (FPO41 + 2# * FPO42 + 3# * FPO43) / 30.973762 / 1000# / 1000# !'eqH+/L per mgP/m^3
                    Pcharge = (FPO41 + 2.0_r64 * FPO42 + 3.0_r64 * FPO43) !avg charge of P species for multiplier for ralkbp

                    SELECT CASE (Rates%IkoxC) !'low O2 inhibition of C oxidation by floating heterotrophs
                      CASE (1)
                        fcarb = c(i, IDX_DISSOLVED_OXYGEN, 2) / (Rates%Ksocf + c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (2)
                        fcarb = 1 - Exp(-Rates%Ksocf * c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (3)
                        fcarb = c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / (Rates%Ksocf + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    END SELECT

                    !gp 08-Feb-06
                    !SELECT CASE (Rates%IkoxCH) !'low O2 inhibition of C oxidation by hyporheic biofilm heterotrophs
                    !CASE (1)
                    ! fcarbH = c(i, IDX_DISSOLVED_OXYGEN, 2) / (Rates%kinhcH + c(i, IDX_DISSOLVED_OXYGEN, 2))
                    !CASE (2)
                    ! fcarbH = 1 - Exp(-Rates%kinhcH * c(i, IDX_DISSOLVED_OXYGEN, 2))
                    !CASE (3)
                    ! fcarbH = c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / (Rates%kinhcH + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    !END SELECT
                    SELECT CASE (Rates%IkoxCH) !'low O2 inhibition of C oxidation by hyporheic biofilm heterotrophs
                      CASE (1)
                        fcarbH = c(i, IDX_DISSOLVED_OXYGEN, 2) / (hydrau%reach(i)%kinhcH + c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (2)
                        fcarbH = 1 - Exp(-hydrau%reach(i)%kinhcH * c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (3)
                        fcarbH = c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / &
                            (hydrau%reach(i)%kinhcH + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxN) !'low O2 inhibition of nitrification
                      CASE (1)
                        fnitr = c(i, IDX_DISSOLVED_OXYGEN, 2) / (Rates%Ksona + c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (2)
                        fnitr = 1 - Exp(-Rates%Ksona * c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (3)
                        fnitr = c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / (Rates%Ksona + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxDN) !'low O2 enhancement of denitrification
                      CASE (1)
                        fdenitr = 1 - c(i, IDX_DISSOLVED_OXYGEN, 2) / (Rates%Ksodn + c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (2)
                        fdenitr = Exp(-Rates%Ksodn * c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (3)
                        fdenitr = 1 - c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / (Rates%Ksodn + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxP) !'low O2 inhib of phytoplankton respiration
                      CASE (1)
                        frespp = c(i, IDX_DISSOLVED_OXYGEN, 2) / (Rates%Ksop + c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (2)
                        frespp = 1 - Exp(-Rates%Ksop * c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (3)
                        frespp = c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / (Rates%Ksop + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    END SELECT

                    !'Phytoplankton (11)
                    PhytoResp = frespp * kreaT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_PHYTOPLANKTON, 2)
                    PhytoDeath = kdeaT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_PHYTOPLANKTON, 2)
                    dc(i, IDX_PHYTOPLANKTON, 2) = dc(i, IDX_PHYTOPLANKTON, 2) - PhytoResp - PhytoDeath

                    !'Enhanced oxidation of fast C in the hyporheic sediment zone
                    !'later versions of Q2K may use c(i, IDX_HYPORHEIC_BIOFILM, 2) to represent biomass of the biofilm
                    !'this current vesion of Q2K does not simulate the biofilm biomass as a state variable
                    !'limitation from fast DOC and dissolved oxygen:

                    !gp 08-Feb-06
                    !phint = 1
                    !IF ((Rates%kscH + c(i, IDX_CBOD_FAST, 2)) <= 0) THEN
                    ! phic = 0
                    !ELSE
                    ! phic = c(i, IDX_CBOD_FAST, 2) / (Rates%kscH + c(i, IDX_CBOD_FAST, 2))
                    !END IF
                    !IF (phic < phint) THEN
                    ! phint = phic
                    !END IF
                    phint = 1
                    IF ((hydrau%reach(i)%kscH + c(i, IDX_CBOD_FAST, 2)) <= 0) THEN
                        phic = 0
                    ELSE
                        phic = c(i, IDX_CBOD_FAST, 2) / (hydrau%reach(i)%kscH + c(i, IDX_CBOD_FAST, 2))
                    END IF
                    IF (phic < phint) THEN
                        phint = phic
                    END IF

                    !'O2 inhibition of heterotroph growth
                    IF (fcarbH < phint) THEN
                        phint = fcarbH
                    END IF
                    !'convert limited zero-order enhanced hyporheic oxygen flux
                    !'or first-order enhanced hyporheic fast CBOD oxidation
                    !'to equivalent units of gD/d of net growth of heterotrophic bacteria biofilm
                    IF (Rates%typeH == "First-order") THEN
                        !'convert limited first-order d^-1 to gD/d of net growth
                        HeteroGrow = (Rates%adc / Rates%roc) * phint * kgaHT(i) * &
                            hydrau%reach(i)%HypoPoreVol * c(i, IDX_CBOD_FAST, 2)
                    ELSE
                        !'convert limited zero-order gO2/m^2/d to gD/d of net growth
                        HeteroGrow = (Rates%adc / Rates%roc) * phint * kgaHT(i) * hydrau%reach(i)%Ast
                    END IF

                    !'detritus (12)
                    DetrDiss = kdtT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_PARTICULATE_ORG_MAT, 2)
                    dc(i, IDX_PARTICULATE_ORG_MAT, 2) = dc(i, IDX_PARTICULATE_ORG_MAT, 2) - DetrDiss
                    dc(i, IDX_PARTICULATE_ORG_MAT, 2) = dc(i, IDX_PARTICULATE_ORG_MAT, 2) + Rates%ada * PhytoDeath

                    !'Organic Nitrogen (6)
                    OrgNHydr = khnT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_ORGANIC_N, 2)
                    dc(i, IDX_ORGANIC_N, 2) = dc(i, IDX_ORGANIC_N, 2) + Rates%ana * PhytoDeath
                    dc(i, IDX_ORGANIC_N, 2) = dc(i, IDX_ORGANIC_N, 2) - OrgNHydr

                    !'Ammonium Nitrogen (7)

                    !gp 03-Dec-09
                    !NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_AMMONIA_N, 2)
                    NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%HypoPoreVol * Fi * c(i, IDX_AMMONIA_N, 2)

                    dc(i, IDX_AMMONIA_N, 2) = dc(i, IDX_AMMONIA_N, 2) + OrgNHydr
                    dc(i, IDX_AMMONIA_N, 2) = dc(i, IDX_AMMONIA_N, 2) - NH4Nitrif
                    dc(i, IDX_AMMONIA_N, 2) = dc(i, IDX_AMMONIA_N, 2) + Rates%ana * PhytoResp

                    !'Nitrate Nitrogen (8)
                    Denitr = c(i, IDX_CBOD_FAST, 2) / (0.1 + c(i, IDX_CBOD_FAST, 2)) * fdenitr * &
                        kiT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_NOX_N, 2)
                    dc(i, IDX_NOX_N, 2) = dc(i, IDX_NOX_N, 2) + NH4Nitrif
                    dc(i, IDX_NOX_N, 2) = dc(i, IDX_NOX_N, 2) - Denitr

                    !'Organic Phosphorus (9)
                    OrgPHydr = khpT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_ORGANIC_P, 2)
                    dc(i, IDX_ORGANIC_P, 2) = dc(i, IDX_ORGANIC_P, 2) + Rates%apa * PhytoDeath
                    dc(i, IDX_ORGANIC_P, 2) = dc(i, IDX_ORGANIC_P, 2) - OrgPHydr

                    !'Inorganic Phosphorus (10)
                    dc(i, IDX_SOLUBLE_REACTIVE_P, 2) = dc(i, IDX_SOLUBLE_REACTIVE_P, 2) + OrgPHydr
                    dc(i, IDX_SOLUBLE_REACTIVE_P, 2) = dc(i, IDX_SOLUBLE_REACTIVE_P, 2) + Rates%apa * PhytoResp

                    !'CBOD Slow (4)
                    CBODsHydr = khcT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_CBOD_SLOW, 2)
                    dc(i, IDX_CBOD_SLOW, 2) = dc(i, IDX_CBOD_SLOW, 2) + Rates%roc * (1 / Rates%adc) * DetrDiss
                    dc(i, IDX_CBOD_SLOW, 2) = dc(i, IDX_CBOD_SLOW, 2) - CBODsHydr
                    CBODsOxid = fcarb * kdcsT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_CBOD_SLOW, 2)
                    dc(i, IDX_CBOD_SLOW, 2) = dc(i, IDX_CBOD_SLOW, 2) - CBODsOxid

                    !'CBOD Fast (5)
                    CBODfOxid = fcarb * kdcT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_CBOD_FAST, 2)
                    CBODfOxid = CBODfOxid + Rates%roc * (1 / Rates%adc) * &
                        HeteroGrow !'enhanced CBOD oxidation from hyporheic biofilm growth
                    dc(i, IDX_CBOD_FAST, 2) = dc(i, IDX_CBOD_FAST, 2) + CBODsHydr
                    dc(i, IDX_CBOD_FAST, 2) = dc(i, IDX_CBOD_FAST, 2) - CBODfOxid
                    dc(i, IDX_CBOD_FAST, 2) = dc(i, IDX_CBOD_FAST, 2) - Rates%rondn * Denitr

                    !gp 08-Dec-04
                    !GENERIC CONSTITUENT or COD (14)
                    dc(i, IDX_GENERIC_CONST, 2) = dc(i, IDX_GENERIC_CONST, 2) - kgenT(i) * &
                        hydrau%reach(i)%HypoPoreVol * c(i, IDX_GENERIC_CONST, 2)
                    IF (Rates%useGenericAsCOD == "Yes") THEN
                        CODoxid = kgenT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_GENERIC_CONST, 2)
                    ELSE
                        CODoxid = 0
                    END IF

                    !'Dissolved Oxygen (3)
                    dc(i, IDX_DISSOLVED_OXYGEN, 2) = dc(i, IDX_DISSOLVED_OXYGEN, 2) - CBODsOxid - CBODfOxid
                    dc(i, IDX_DISSOLVED_OXYGEN, 2) = dc(i, IDX_DISSOLVED_OXYGEN, 2) - Rates%ron * NH4Nitrif
                    dc(i, IDX_DISSOLVED_OXYGEN, 2) = dc(i, IDX_DISSOLVED_OXYGEN, 2) - Rates%roa * PhytoResp
                    dc(i, IDX_DISSOLVED_OXYGEN, 2) = dc(i, IDX_DISSOLVED_OXYGEN, 2) - CODoxid !gp 08-Dec-04

                    !'Alkalinity (nv - 2)

                    If (sys%simAlk == "Yes") Then !gp 26-Oct-07

                        !gp 03-Dec-09
                        !dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkaa * PhytoResp * 50000
                        !dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) - Rates%ralkn * NH4Nitrif * 50000
                        !dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkden * Denitr * 50000
                        !dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkbn * OrgNHydr * 50000 + Rates%ralkbp * OrgPHydr * 50000
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkbn * Fi * &
                            Rates%ana * PhytoResp * 50043.45_r64 !'phyto resp of N
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) - Rates%ralkbp * Pcharge * &
                            Rates%apa * PhytoResp * 50043.45_r64 !'phyto resp of P
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) - Rates%ralkbn * 2.0_r64 * &
                            NH4Nitrif * 50043.45_r64 !'nitrification ammonia loss (Fi) and nitrate gain (+1)
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkbn * &
                            Denitr * 50043.45_r64 !'denitrification
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkbn * Fi * &
                            OrgNHydr * 50043.45_r64 !'organic N hydrolysis
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) - Rates%ralkbp * Pcharge * &
                            OrgPHydr * 50043.45_r64 !'organic P hydrolysis

                    end if !gp 26-Oct-07

                    !'Inorganic Carbon (nv - 1)
                    dc(i, IDX_TOTAL_INORGANIC_C, 2) = dc(i, IDX_TOTAL_INORGANIC_C, 2) + Rates%rcca * PhytoResp
                    dc(i, IDX_TOTAL_INORGANIC_C, 2) = dc(i, IDX_TOTAL_INORGANIC_C, 2) + Rates%rcco * CBODfOxid
                    dc(i, IDX_TOTAL_INORGANIC_C, 2) = dc(i, IDX_TOTAL_INORGANIC_C, 2) + Rates%rcco * CBODsOxid

                    !'PATHOGEN
                    dc(i, IDX_PATHOGEN, 2) = dc(i, IDX_PATHOGEN, 2) - kpathT(i) * &
                        hydrau%reach(i)%HypoPoreVol * c(i, IDX_PATHOGEN, 2)

                END DO

                !'
                !' ------------------------------------------------------------------------
                !' --- convert hyporheic pore water deriviatives to concentration units ---
                !' ------------------------------------------------------------------------
                !'

                DO i = 1, nr
                    DO k = 1, nv - 1
                        dc(i, k, 2) = dc(i, k, 2) / hydrau%reach(i)%HypoPoreVol
                    END DO
                END DO


                !gp 15-Nov-04
                !
                ! -----------------------------------------------------
                ! -----------------------------------------------------
                ! --- Level 2 growth kinetics for hyporheic biofilm ---
                ! -----------------------------------------------------
                ! -----------------------------------------------------
                !

              CASE ('Level 2') !gp 15-Nov-04

                !
                ! ------------------------
                ! --- calc and save pH ---
                ! ------------------------
                !

                DO i=0, nr
                    call ph_solver(sys%imethph, ph, c(i, IDX_TOTAL_INORGANIC_C, 2), Te(i, 2), &
                        c(i, IDX_ALKALINITY, 2), c(i, IDX_CONDUCTIVITY, 2))
                    pHs(i, 2) = pH
                    CALL ChemRates(Te(i, 2), K1, K2, KW, Kh, c(i, IDX_CONDUCTIVITY, 2))
                    K1s(i, 2) = K1; K2s(i, 2) = K2; Khs(i, 2) = Kh
                END DO

                !
                ! -----------------------------------------------------
                ! --- temperature adjustment of rates and constants ---
                ! -----------------------------------------------------
                !

                DO i = 1, nr

                    !gp 08-Feb-06
                    !khcT(i) = Rates%khc * Rates%tkhc ** (Te(i, 2) - 20)
                    !kdcsT(i) = Rates%kdcs * Rates%tkdcs ** (Te(i, 2) - 20)
                    !kdcT(i) = Rates%kdc * Rates%tkdc ** (Te(i, 2) - 20)
                    !khnT(i) = Rates%khn * Rates%tkhn ** (Te(i, 2) - 20)
                    !khpT(i) = Rates%khp * Rates%tkhp ** (Te(i, 2) - 20)
                    !knT(i) = Rates%kn * Rates%tkn ** (Te(i, 2) - 20)
                    !kiT(i) = Rates%ki * Rates%tki ** (Te(i, 2) - 20)
                    !os(i, 2)=oxsat(Te(i, 2), hydrau%reach(i)%elev)
                    !kdeaT(i) = Rates%kdea * Rates%tkdea ** (Te(i, 2) - 20)
                    !kreaT(i) = Rates%krea * Rates%tkrea ** (Te(i, 2) - 20)
                    !kdtT(i) = Rates%kdt * Rates%tkdt ** (Te(i, 2) - 20)
                    !kpathT(i) = Rates%kpath * Rates%tkpath ** (Te(i, 2) - 20)
                    !kgaHT(i) = Rates%kgaH * Rates%tkgaH ** (Te(i, 2) - 20) !enhanced CBOD oxidation rate from hyporheic biofilm
                    !kgaHT(i) = Rates%kgaH * Rates%tkgaH ** (Te(i, 2) - 20) !gp 15-Nov-04 hyporheic biofilm growth
                    !kreaHT(i) = Rates%kreaH * Rates%tkreaH ** (Te(i, 2) - 20) !gp 15-Nov-04 hyporheic biofilm respiration
                    !kdeaHT(i) = Rates%kdeaH * Rates%tkdeaH ** (Te(i, 2) - 20) !gp 15-Nov-04 hyporheic biofilm death
                    !kgenT(i) = Rates%kgen * Rates%tkgen ** (Te(i, 2) - 20) !gp 30-Nov-04
                    khcT(i) = hydrau%reach(i)%khc * Rates%tkhc ** (Te(i, 2) - 20)
                    kdcsT(i) = hydrau%reach(i)%kdcs * Rates%tkdcs ** (Te(i, 2) - 20)
                    kdcT(i) = hydrau%reach(i)%kdc * Rates%tkdc ** (Te(i, 2) - 20)
                    khnT(i) = hydrau%reach(i)%khn * Rates%tkhn ** (Te(i, 2) - 20)
                    khpT(i) = hydrau%reach(i)%khp * Rates%tkhp ** (Te(i, 2) - 20)
                    knT(i) = hydrau%reach(i)%kn * Rates%tkn ** (Te(i, 2) - 20)
                    kiT(i) = hydrau%reach(i)%ki * Rates%tki ** (Te(i, 2) - 20)
                    os(i, 2)=oxygen_saturation(Te(i, 2), hydrau%reach(i)%elev)
                    kdeaT(i) = hydrau%reach(i)%kdea * Rates%tkdea ** (Te(i, 2) - 20)
                    kreaT(i) = hydrau%reach(i)%krea * Rates%tkrea ** (Te(i, 2) - 20)
                    kdtT(i) = hydrau%reach(i)%kdt * Rates%tkdt ** (Te(i, 2) - 20)
                    kpathT(i) = hydrau%reach(i)%kpath * Rates%tkpath ** (Te(i, 2) - 20)
                    kgaHT(i) = hydrau%reach(i)%kgaH * Rates%tkgaH ** (Te(i, 2) - 20)
                    kreaHT(i) = hydrau%reach(i)%kreaH * Rates%tkreaH ** (Te(i, 2) - 20)
                    kdeaHT(i) = hydrau%reach(i)%kdeaH * Rates%tkdeaH ** (Te(i, 2) - 20)
                    kgenT(i) = hydrau%reach(i)%kgen * Rates%tkgen ** (Te(i, 2) - 20)

                END DO

                !
                ! ---------------------------------------------------------------------------------------------
                ! --- mass transfer between water column and hyporheic zone by bulk hyporheic exchange flow ---
                ! ---------------------------------------------------------------------------------------------
                !

                DO i = 1, nr
                    DO k = 1, nv-1
                        !gp un-comment next two lines to activate mass transfer between water column and hyporheic pore water
                        dc(i, k, 1) = dc(i, k, 1) + hydrau%reach(i)%EhyporheicCMD * (c(i, k, 2) - c(i, k, 1))
                        dc(i, k, 2) = hydrau%reach(i)%EhyporheicCMD * (c(i, k, 1) - c(i, k, 2))
                    END DO
                    !store hyporheic fluxes for output (positive flux is source to the water column from the hyporheic zone)
                    HypoFluxDO(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_DISSOLVED_OXYGEN, 2) - c(i, IDX_DISSOLVED_OXYGEN, 1)) / hydrau%reach(i)%Ast
                    ! DO gO2/m^2/d
                    HypoFluxCBOD(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_CBOD_FAST, 2) - c(i, IDX_CBOD_FAST, 1)) / hydrau%reach(i)%Ast
                    ! fast CBOD gO2/m^2/d
                    HypoFluxNH4(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_AMMONIA_N, 2) - c(i, IDX_AMMONIA_N, 1)) / hydrau%reach(i)%Ast
                    ! NH4 mgN/m^2/d
                    HypoFluxNO3(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_NOX_N, 2) - c(i, IDX_NOX_N, 1)) / hydrau%reach(i)%Ast
                    ! NO3 mgN/m^2/d
                    HypoFluxSRP(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_SOLUBLE_REACTIVE_P, 2) - c(i, IDX_SOLUBLE_REACTIVE_P, 1)) / hydrau%reach(i)%Ast
                    ! SRP mgP/m^2/d
                    HypoFluxIC(i) = hydrau%reach(i)%EhyporheicCMD * &
                        (c(i, IDX_TOTAL_INORGANIC_C, 2) - c(i, IDX_TOTAL_INORGANIC_C, 1)) / &
                        hydrau%reach(i)%Ast / Rates%rccc !cT gC/m^2/d
                END DO

                !
                ! -----------------------------------------------------
                ! --- mass derivatives for the hyporheic pore water ---
                ! -----------------------------------------------------
                !

                DO i = 1, nr

                    ! --- pH dependent varialbes for later use in deriv calcs in this i loop

                    !'gp 03-Dec-09
                    !'fraction of ionized ammonia
                    hh = 10.0_r64 ** (-pHs(i, 2))
                    Kamm = 10.0_r64 ** (-(0.09018_r64 + 2729.92_r64 / (Te(i, 2) + 273.15_r64)))
                    Fi = hh / (hh + Kamm)
                    !'phosphate fractions of H2PO4- (FPO41), HPO4-- (FPO42), and PO4--- (FPO43)
                    KPO41 = 10.0_r64 ** (-2.15_r64)
                    KPO42 = 10.0_r64 ** (-7.2_r64)
                    KPO43 = 10.0_r64 ** (-12.35_r64)
                    DPO = 1.0_r64 / (hh ** 3.0_r64 + KPO41 * hh ** 2.0_r64 + KPO41 * KPO42 * hh + KPO41 * KPO42 * KPO43)
                    FPO41 = KPO41 * hh ** 2.0_r64 * DPO !'fraction of phosphate as H2PO4-
                    FPO42 = KPO41 * KPO42 * hh * DPO !'fraction of phosphate as HPO4--
                    FPO43 = KPO41 * KPO42 * KPO43 * DPO !'fraction of phosphate as PO4---
                    !'stoichiometric conversion factors (some are pH dependent)
                    !ralkbn = 1# / 14.0067 / 1000# / 1000# !'eqH+/L per mgN/m^3
                    !ralkbp = (FPO41 + 2# * FPO42 + 3# * FPO43) / 30.973762 / 1000# / 1000# !'eqH+/L per mgP/m^3
                    Pcharge = (FPO41 + 2.0_r64 * FPO42 + 3.0_r64 * FPO43) !avg charge of P species for multiplier for ralkbp

                    SELECT CASE (Rates%IkoxC) !'low O2 inhibition of C oxidation by floating heterotrophs
                      CASE (1)
                        fcarb = c(i, IDX_DISSOLVED_OXYGEN, 2) / (Rates%Ksocf + c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (2)
                        fcarb = 1 - Exp(-Rates%Ksocf * c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (3)
                        fcarb = c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / (Rates%Ksocf + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    END SELECT

                    !gp 08-Feb-06
                    !SELECT CASE (Rates%IkoxCH) !'low O2 inhibition of C oxidation by hyporheic biofilm heterotrophs
                    !CASE (1)
                    ! fcarbH = c(i, IDX_DISSOLVED_OXYGEN, 2) / (Rates%kinhcH + c(i, IDX_DISSOLVED_OXYGEN, 2))
                    !CASE (2)
                    ! fcarbH = 1 - Exp(-Rates%kinhcH * c(i, IDX_DISSOLVED_OXYGEN, 2))
                    !CASE (3)
                    ! fcarbH = c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / (Rates%kinhcH + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    !END SELECT
                    SELECT CASE (Rates%IkoxCH) !'low O2 inhibition of C oxidation by hyporheic biofilm heterotrophs
                      CASE (1)
                        fcarbH = c(i, IDX_DISSOLVED_OXYGEN, 2) / (hydrau%reach(i)%kinhcH + c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (2)
                        fcarbH = 1 - Exp(-hydrau%reach(i)%kinhcH * c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (3)
                        fcarbH = c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / &
                            (hydrau%reach(i)%kinhcH + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxN) !'low O2 inhibition of nitrification
                      CASE (1)
                        fnitr = c(i, IDX_DISSOLVED_OXYGEN, 2) / (Rates%Ksona + c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (2)
                        fnitr = 1 - Exp(-Rates%Ksona * c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (3)
                        fnitr = c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / (Rates%Ksona + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxDN) !'low O2 enhancement of denitrification
                      CASE (1)
                        fdenitr = 1 - c(i, IDX_DISSOLVED_OXYGEN, 2) / (Rates%Ksodn + c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (2)
                        fdenitr = Exp(-Rates%Ksodn * c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (3)
                        fdenitr = 1 - c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / (Rates%Ksodn + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    END SELECT

                    SELECT CASE (Rates%IkoxP) !'low O2 inhib of phytoplankton respiration
                      CASE (1)
                        frespp = c(i, IDX_DISSOLVED_OXYGEN, 2) / (Rates%Ksop + c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (2)
                        frespp = 1 - Exp(-Rates%Ksop * c(i, IDX_DISSOLVED_OXYGEN, 2))
                      CASE (3)
                        frespp = c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2 / (Rates%Ksop + c(i, IDX_DISSOLVED_OXYGEN, 2) ** 2)
                    END SELECT

                    !'Phytoplankton (11)
                    PhytoResp = frespp * kreaT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_PHYTOPLANKTON, 2)
                    PhytoDeath = kdeaT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_PHYTOPLANKTON, 2)
                    dc(i, IDX_PHYTOPLANKTON, 2) = dc(i, IDX_PHYTOPLANKTON, 2) - PhytoResp - PhytoDeath

                    !gp 15-Nov-04
                    !Attached heterotrophic bacteria (nv)
                    !nutrient limitation

                    !gp 03-Dec-09
                    !If (hydrau%reach(i)%ksnH + c(i, IDX_AMMONIA_N, 2) + c(i, IDX_NOX_N, 2) <= 0) Then
                    ! phint = 0
                    !Else
                    ! phint = (c(i, IDX_AMMONIA_N, 2) + c(i, IDX_NOX_N, 2)) / (hydrau%reach(i)%ksnH + c(i, IDX_AMMONIA_N, 2) + c(i, IDX_NOX_N, 2))
                    !End If
                    !If (hydrau%reach(i)%kspH + c(i, IDX_SOLUBLE_REACTIVE_P, 2) <= 0) Then
                    ! phip = 0
                    !Else
                    ! phip = c(i, IDX_SOLUBLE_REACTIVE_P, 2) / (hydrau%reach(i)%kspH + c(i, IDX_SOLUBLE_REACTIVE_P, 2))
                    !End If
                    !If (phip < phint) phint = phip
                    If (hydrau%reach(i)%ksnH + Fi * c(i, IDX_AMMONIA_N, 2) + c(i, IDX_NOX_N, 2) <= 0) Then
                        phint = 0
                    Else
                        phint = (Fi * c(i, IDX_AMMONIA_N, 2) + c(i, IDX_NOX_N, 2)) / &
                            (hydrau%reach(i)%ksnH + Fi * c(i, IDX_AMMONIA_N, 2) + c(i, IDX_NOX_N, 2))
                    End If
                    If (hydrau%reach(i)%kspH + c(i, IDX_SOLUBLE_REACTIVE_P, 2) <= 0) Then
                        phip = 0
                    Else
                        phip = c(i, IDX_SOLUBLE_REACTIVE_P, 2) / (hydrau%reach(i)%kspH + c(i, IDX_SOLUBLE_REACTIVE_P, 2))
                    End If
                    If (phip < phint) phint = phip

                    !limitation of enhanced oxidation of fast DOC from fast DOC

                    If (hydrau%reach(i)%kscH + c(i, IDX_CBOD_FAST, 2) <= 0) Then
                        phic = 0
                    Else
                        phic = c(i, IDX_CBOD_FAST, 2) / (hydrau%reach(i)%kscH + c(i, IDX_CBOD_FAST, 2))
                    End If

                    If (phic < phint) phint = phic
                    !O2 inhibition of enhanced oxidation of fast DOC
                    If (fcarbH < phint) phint = fcarbH
                    !optional first-order growth kinetics
                    If (Rates%typeH == "First-order") Then

                        !gp 08-Feb-06
                        !HeteroGrow = phint * kgaHT(i) * hydrau%reach(i)%Ast * c(i, IDX_HYPORHEIC_BIOFILM, 2) * (1 - c(i, IDX_HYPORHEIC_BIOFILM, 2) / Rates%ahmax) !use first-order growth kinetics
                        HeteroGrow = phint * kgaHT(i) * hydrau%reach(i)%Ast * c(i, IDX_HYPORHEIC_BIOFILM, 2) &
                            * (1 - c(i, IDX_HYPORHEIC_BIOFILM, 2) / hydrau%reach(i)%ahmax) !use first-order growth kinetics

                    Else
                        HeteroGrow = (Rates%adc / Rates%roc) * phint * kgaHT(i) * &
                            hydrau%reach(i)%Ast !use zero-order, kgaH is gO2/m2/d
                    End If
                    HeteroResp = fcarbH * kreaHT(i) * hydrau%reach(i)%Ast * &
                        c(i, IDX_HYPORHEIC_BIOFILM, 2) !first-order respiration with O2 limitation
                    HeteroDeath = kdeaHT(i) * hydrau%reach(i)%Ast * c(i, IDX_HYPORHEIC_BIOFILM, 2) !first-order death
                    dc(i, IDX_HYPORHEIC_BIOFILM, 2) = HeteroGrow - HeteroResp - HeteroDeath !units gD/d deriv for biofilm biomass

                    !ammonium preference

                    !gp 03-Dec-09
                    !prefamH = 0
                    !If (c(i, IDX_AMMONIA_N, 2) + c(i, IDX_NOX_N, 2) > 0) Then
                    ! prefamH = c(i, IDX_AMMONIA_N, 2) * c(i, IDX_NOX_N, 2) / (hydrau%reach(i)%khnxH + c(i, IDX_AMMONIA_N, 2)) / (hydrau%reach(i)%khnxH + c(i, IDX_NOX_N, 2)) &
                    ! + c(i, IDX_AMMONIA_N, 2) * hydrau%reach(i)%khnxH / (c(i, IDX_AMMONIA_N, 2) + c(i, IDX_NOX_N, 2)) / (hydrau%reach(i)%khnxH + c(i, IDX_NOX_N, 2))
                    !End If
                    prefamH = 0
                    If (Fi * c(i, IDX_AMMONIA_N, 2) + c(i, IDX_NOX_N, 2) > 0) Then
                        prefamH = Fi * c(i, IDX_AMMONIA_N, 2) * c(i, IDX_NOX_N, 2) / &
                            (hydrau%reach(i)%khnxH + Fi * c(i, IDX_AMMONIA_N, 2)) / &
                            (hydrau%reach(i)%khnxH + c(i, IDX_NOX_N, 2)) &
                            + Fi * c(i, IDX_AMMONIA_N, 2) * hydrau%reach(i)%khnxH / &
                            (Fi * c(i, IDX_AMMONIA_N, 2) + c(i, IDX_NOX_N, 2)) / &
                            (hydrau%reach(i)%khnxH + c(i, IDX_NOX_N, 2))
                    End If

                    !'detritus (12)
                    DetrDiss = kdtT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_PARTICULATE_ORG_MAT, 2)
                    dc(i, IDX_PARTICULATE_ORG_MAT, 2) = dc(i, IDX_PARTICULATE_ORG_MAT, 2) - DetrDiss
                    dc(i, IDX_PARTICULATE_ORG_MAT, 2) = dc(i, IDX_PARTICULATE_ORG_MAT, 2) + Rates%ada * PhytoDeath
                    dc(i, IDX_PARTICULATE_ORG_MAT, 2) = dc(i, IDX_PARTICULATE_ORG_MAT, 2) + HeteroDeath !gp 22-Nov-04

                    !'Organic Nitrogen (6)
                    OrgNHydr = khnT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_ORGANIC_N, 2)
                    dc(i, IDX_ORGANIC_N, 2) = dc(i, IDX_ORGANIC_N, 2) + Rates%ana * PhytoDeath
                    dc(i, IDX_ORGANIC_N, 2) = dc(i, IDX_ORGANIC_N, 2) - OrgNHydr
                    dc(i, IDX_ORGANIC_N, 2) = dc(i, IDX_ORGANIC_N, 2) + HeteroDeath * Rates%ana / Rates%ada !gp 15-Nov-04

                    !'Ammonium Nitrogen (7)

                    !gp 03-Dec-09
                    !NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_AMMONIA_N, 2)
                    NH4Nitrif = fnitr * knT(i) * hydrau%reach(i)%HypoPoreVol * Fi * c(i, IDX_AMMONIA_N, 2)

                    dc(i, IDX_AMMONIA_N, 2) = dc(i, IDX_AMMONIA_N, 2) + OrgNHydr
                    dc(i, IDX_AMMONIA_N, 2) = dc(i, IDX_AMMONIA_N, 2) - NH4Nitrif
                    dc(i, IDX_AMMONIA_N, 2) = dc(i, IDX_AMMONIA_N, 2) + Rates%ana * PhytoResp
                    dc(i, IDX_AMMONIA_N, 2) = dc(i, IDX_AMMONIA_N, 2) + Rates%anc / Rates%adc * HeteroResp !gp 15-Nov-04
                    dc(i, IDX_AMMONIA_N, 2) = dc(i, IDX_AMMONIA_N, 2) - Rates%anc / Rates%adc * prefamH * HeteroGrow !gp 15-Nov-04

                    !'Nitrate Nitrogen (8)
                    Denitr = c(i, IDX_CBOD_FAST, 2) / (0.1 + c(i, IDX_CBOD_FAST, 2)) * fdenitr * &
                        kiT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_NOX_N, 2)
                    dc(i, IDX_NOX_N, 2) = dc(i, IDX_NOX_N, 2) + NH4Nitrif
                    dc(i, IDX_NOX_N, 2) = dc(i, IDX_NOX_N, 2) - Denitr
                    dc(i, IDX_NOX_N, 2) = dc(i, IDX_NOX_N, 2) - Rates%anc / Rates%adc * (1 - prefamH) * HeteroGrow !gp 15-Nov-04

                    !'Organic Phosphorus (9)
                    OrgPHydr = khpT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_ORGANIC_P, 2)
                    dc(i, IDX_ORGANIC_P, 2) = dc(i, IDX_ORGANIC_P, 2) + Rates%apa * PhytoDeath
                    dc(i, IDX_ORGANIC_P, 2) = dc(i, IDX_ORGANIC_P, 2) - OrgPHydr
                    dc(i, IDX_ORGANIC_P, 2) = dc(i, IDX_ORGANIC_P, 2) + HeteroDeath * Rates%apa / Rates%ada !gp 15-Nov-04

                    !'Inorganic Phosphorus (10)
                    dc(i, IDX_SOLUBLE_REACTIVE_P, 2) = dc(i, IDX_SOLUBLE_REACTIVE_P, 2) + OrgPHydr
                    dc(i, IDX_SOLUBLE_REACTIVE_P, 2) = dc(i, IDX_SOLUBLE_REACTIVE_P, 2) + Rates%apa * PhytoResp
                    dc(i, IDX_SOLUBLE_REACTIVE_P, 2) = dc(i, IDX_SOLUBLE_REACTIVE_P, 2) + Rates%apc / Rates%adc * &
                        HeteroResp !gp 15-Nov-04
                    dc(i, IDX_SOLUBLE_REACTIVE_P, 2) = dc(i, IDX_SOLUBLE_REACTIVE_P, 2) - Rates%apc / Rates%adc * &
                        HeteroGrow !gp 15-Nov-04

                    !'CBOD Slow (4)
                    CBODsHydr = khcT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_CBOD_SLOW, 2)
                    dc(i, IDX_CBOD_SLOW, 2) = dc(i, IDX_CBOD_SLOW, 2) + Rates%roc * (1 / Rates%adc) * DetrDiss
                    dc(i, IDX_CBOD_SLOW, 2) = dc(i, IDX_CBOD_SLOW, 2) - CBODsHydr
                    CBODsOxid = fcarb * kdcsT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_CBOD_SLOW, 2)
                    dc(i, IDX_CBOD_SLOW, 2) = dc(i, IDX_CBOD_SLOW, 2) - CBODsOxid

                    !'CBOD Fast (5)
                    CBODfOxid = fcarb * kdcT(i) * hydrau%reach(i)%HypoPoreVol * &
                        c(i, IDX_CBOD_FAST, 2) !oxidation by suspended bacteria
                    CBODfOxid = CBODfOxid + Rates%roc * (1 / Rates%adc) * &
                        HeteroGrow !enhanced CBOD oxidation from hyporheic biofilm growth
                    dc(i, IDX_CBOD_FAST, 2) = dc(i, IDX_CBOD_FAST, 2) + CBODsHydr
                    dc(i, IDX_CBOD_FAST, 2) = dc(i, IDX_CBOD_FAST, 2) - CBODfOxid
                    dc(i, IDX_CBOD_FAST, 2) = dc(i, IDX_CBOD_FAST, 2) - Rates%rondn * Denitr

                    !gp 08-Dec-04
                    !GENERIC CONSTITUENT or COD (14)
                    dc(i, IDX_GENERIC_CONST, 2) = dc(i, IDX_GENERIC_CONST, 2) - kgenT(i) * &
                        hydrau%reach(i)%HypoPoreVol * c(i, IDX_GENERIC_CONST, 2)
                    IF (Rates%useGenericAsCOD == "Yes") THEN
                        CODoxid = kgenT(i) * hydrau%reach(i)%HypoPoreVol * c(i, IDX_GENERIC_CONST, 2)
                    ELSE
                        CODoxid = 0
                    END IF

                    !'Dissolved Oxygen (3)
                    dc(i, IDX_DISSOLVED_OXYGEN, 2) = dc(i, IDX_DISSOLVED_OXYGEN, 2) - CBODsOxid - &
                        CBODfOxid !CBODfOxid includes oxidation by hyporheic biofilm
                    dc(i, IDX_DISSOLVED_OXYGEN, 2) = dc(i, IDX_DISSOLVED_OXYGEN, 2) - Rates%ron * NH4Nitrif
                    dc(i, IDX_DISSOLVED_OXYGEN, 2) = dc(i, IDX_DISSOLVED_OXYGEN, 2) - Rates%roa * PhytoResp
                    dc(i, IDX_DISSOLVED_OXYGEN, 2) = dc(i, IDX_DISSOLVED_OXYGEN, 2) - Rates%roc / Rates%adc * &
                        HeteroResp !gp 15-Nov-04
                    dc(i, IDX_DISSOLVED_OXYGEN, 2) = dc(i, IDX_DISSOLVED_OXYGEN, 2) - CODoxid !gp 08-Dec-04

                    !'Alkalinity (nv - 2)

                    If (sys%simAlk == "Yes") Then !gp 26-Oct-07

                        !gp 03-Dec-09
                        !dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkaa * PhytoResp * 50000
                        !dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) - Rates%ralkn * NH4Nitrif * 50000
                        !dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkden * Denitr * 50000
                        !dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkbn * OrgNHydr * 50000 + Rates%ralkbp * OrgPHydr * 50000
                        !dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) - Rates%ralkda * HeteroGrow * prefamH * 50000 !gp 15-Nov-04
                        !dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkdn * HeteroGrow * (1 - prefamH) * 50000 !gp 15-Nov-04
                        !dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkda * HeteroResp * 50000 !gp 15-Nov-04
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkbn * Fi * &
                            Rates%ana * PhytoResp * 50043.45_r64 !'phyto resp of N
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) - Rates%ralkbp * Pcharge * &
                            Rates%apa * PhytoResp * 50043.45_r64 !'phyto resp of P
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) - Rates%ralkbn * 2.0_r64 * &
                            NH4Nitrif * 50043.45_r64 !'nitrification ammonia loss (Fi) and nitrate gain (+1)
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkbn * &
                            Denitr * 50043.45_r64 !'denitrification
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkbn * Fi * &
                            OrgNHydr * 50043.45_r64 !'organic N hydrolysis
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) - Rates%ralkbp * Pcharge * &
                            OrgPHydr * 50043.45_r64 !'organic P hydrolysis
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) - Rates%ralkbn * Rates%anc / Rates%adc * &
                            prefamH * HeteroGrow * 50043.45_r64 !'hetero uptake of ammonia
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkbn * Rates%anc / Rates%adc * &
                            (1 - prefamH) * HeteroGrow * 50043.45_r64 !'hetero uptake of nitrate
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkbp * Pcharge * Rates%apc / &
                            Rates%adc * HeteroGrow * 50043.45_r64 !'hetero uptake of SRP
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) + Rates%ralkbn * Fi * Rates%anc / &
                            Rates%adc * HeteroResp * 50043.45_r64 !'hetero resp of ammonia
                        dc(i, IDX_ALKALINITY, 2) = dc(i, IDX_ALKALINITY, 2) - Rates%ralkbp * Pcharge * Rates%apc / &
                            Rates%adc * HeteroResp * 50043.45_r64 !'hetero resp of SRP

                    end if !gp 26-Oct-07

                    !'Inorganic Carbon (nv - 1)
                    dc(i, IDX_TOTAL_INORGANIC_C, 2) = dc(i, IDX_TOTAL_INORGANIC_C, 2) + Rates%rcca * PhytoResp
                    dc(i, IDX_TOTAL_INORGANIC_C, 2) = dc(i, IDX_TOTAL_INORGANIC_C, 2) + Rates%rcco * CBODfOxid
                    dc(i, IDX_TOTAL_INORGANIC_C, 2) = dc(i, IDX_TOTAL_INORGANIC_C, 2) + Rates%rcco * CBODsOxid
                    dc(i, IDX_TOTAL_INORGANIC_C, 2) = dc(i, IDX_TOTAL_INORGANIC_C, 2) + Rates%rccd * HeteroResp !gp 15-Nov-04

                    !'PATHOGEN
                    dc(i, IDX_PATHOGEN, 2) = dc(i, IDX_PATHOGEN, 2) - kpathT(i) * &
                        hydrau%reach(i)%HypoPoreVol * c(i, IDX_PATHOGEN, 2)

                END DO

                !'
                !' ------------------------------------------------------------------------
                !' --- convert hyporheic pore water deriviatives to concentration units ---
                !' ------------------------------------------------------------------------
                !'

                DO i = 1, nr
                    DO k = 1, nv - 1
                        dc(i, k, 2) = dc(i, k, 2) / hydrau%reach(i)%HypoPoreVol
                    END DO
                    dc(i, IDX_HYPORHEIC_BIOFILM, 2) = dc(i, IDX_HYPORHEIC_BIOFILM, 2) / hydrau%reach(i)%Ast !gp 15-Nov-04
                END DO


            END SELECT !'gp end of calculation of hyporheic water quality derivatives


            !convert water column derivatives from mass to concentration units
            call conv_wat_col_deriv_from_mass_to_conc(nr, nv, hydrau, dc, dinb, dipb)

            !gp 03-Feb-05 end of bypass derivs for water quality variables
            !unless 'All' state variables are being simulated

        End If

    END SUBROUTINE

    subroutine conv_wat_col_deriv_from_mass_to_conc(nr_, nv_, hydrau_, dc_, dinb_, dipb_)
        ! convert water column derivatives from mass to concentration units
        integer(i32), intent(in) :: nr_, nv_
        real(r64), intent(inout):: dc_(:, :, :), dinb_(:), dipb_(:)
        type(riverhydraulics_type), intent(in) :: hydrau_

        integer(i32) :: i, k

        do i = 1, nr_
            do k = 1, nv_ - 1
                dc_(i, k, 1) = dc_(i, k, 1) / hydrau_%reach(i)%vol
            end do
            dc_(i, nv_, 1) = dc_(i, nv_, 1) / hydrau_%reach(i)%asb
            dinb_(i) = dinb_(i) / hydrau_%reach(i)%asb
            dipb_(i) = dipb_(i) / hydrau_%reach(i)%asb
        end do
    end Subroutine conv_wat_col_deriv_from_mass_to_conc

end module m_derivs
