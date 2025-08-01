cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 14)

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  if(DEFINED ENV{VITASDK})
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VITASDK}/share/vita.toolchain.cmake" CACHE PATH "toolchain file")
  else()
    message(FATAL_ERROR "Please define VITASDK to point to your SDK path!")
  endif()
endif()

if (CMAKE_SYSTEM_NAME MATCHES "Emscripten")
    include_directories(src/lib/SDL2/include)
endif ()

project(jedi_academy)

# Local files import directories for easier #include
include_directories(./src/code)
include_directories(./src/code/game)
include_directories(./src/code/rd-common)
include_directories(./src/code/rd-vanilla)
include_directories(./src/shared)
include_directories(./src/lib/gsl-lite/include)
include_directories(./src/lib/minizip/include/minizip)
include_directories(./src/lib/minizip/include)

add_compile_definitions(ARCH_STRING="x86_64")
add_compile_definitions(SP_GAME)
add_compile_definitions(_JK2EXE)

if (CMAKE_SYSTEM_NAME MATCHES "Emscripten")
    add_compile_definitions(__linux__)
endif ()

add_executable(
        ${PROJECT_NAME}.elf
        ./src/shared/sdl/sdl_icon.h
        ./src/shared/sdl/sdl_sound.h
        ./src/shared/sdl/sdl_qgl.h
        ./src/shared/sys/sys_public.h
        ./src/shared/sys/sys_local.h
        ./src/shared/sys/con_local.h
        ./src/shared/sys/unix_local.h
        ./src/shared/sys/sys_loadlib.h
        ./src/shared/qcommon/q_platform.h
        ./src/shared/qcommon/q_color.h
        ./src/shared/qcommon/q_string.h
        ./src/shared/qcommon/q_math.h
        ./src/shared/qcommon/safe/gsl.h
        ./src/shared/qcommon/safe/files.h
        ./src/shared/qcommon/safe/string.h
        ./src/shared/qcommon/safe/limited_vector.h
        ./src/shared/qcommon/safe/sscanf.h
        ./src/code/game/bg_public.h
        ./src/code/game/genericparser2.h
        ./src/code/game/hitlocs.h
        ./src/code/game/g_navigator.h
        ./src/code/game/bset.h
        ./src/code/game/g_roff.h
        ./src/code/game/g_nav.h
        ./src/code/game/w_local.h
        ./src/code/game/say.h
        ./src/code/game/g_public.h
        ./src/code/game/teams.h
        ./src/code/game/g_shared.h
        ./src/code/game/statindex.h
        ./src/code/game/events.h
        ./src/code/game/npc_headers.h
        ./src/code/game/ghoul2_shared.h
        ./src/code/game/g_vehicles.h
        ./src/code/game/wp_saber.h
        ./src/code/game/surfaceflags.h
        ./src/code/game/g_local.h
        ./src/code/game/channels.h
        ./src/code/game/weapons.h
        ./src/code/game/dmstates.h
        ./src/code/game/objectives.h
        ./src/code/game/b_local.h
        ./src/code/game/anims.h
        ./src/code/game/Q3_Interface.h
        ./src/code/game/b_public.h
        ./src/code/game/g_functions.h
        ./src/code/game/bstate.h
        ./src/code/game/bg_local.h
        ./src/code/game/g_items.h
        ./src/code/game/fields.h
        ./src/code/game/ai.h
        ./src/code/game/common_headers.h
        ./src/code/rd-vanilla/qgl.h
        ./src/code/rd-vanilla/tr_WorldEffects.h
        ./src/code/rd-vanilla/tr_stl.h
        ./src/code/rd-vanilla/tr_local.h
        ./src/code/rd-vanilla/glext.h
        ./src/code/rd-vanilla/tr_quicksprite.h
        ./src/code/win32/AutoVersion.h
        ./src/code/win32/resource.h
        ./src/code/win32/afxres.h
        ./src/code/Ravl/CBounds.h
        ./src/code/Ravl/CVec.h
        ./src/code/Ravl/CMatrix.h
        ./src/code/ghoul2/G2.h
        ./src/code/ghoul2/ghoul2_gore.h
        ./src/code/cgame/FxUtil.h
        ./src/code/cgame/FxSystem.h
        ./src/code/cgame/cg_media.h
        ./src/code/cgame/strip_objectives.h
        ./src/code/cgame/animtable.h
        ./src/code/cgame/FxPrimitives.h
        ./src/code/cgame/cg_public.h
        ./src/code/cgame/cg_local.h
        ./src/code/cgame/FxScheduler.h
        ./src/code/cgame/cg_headers.h
        ./src/code/cgame/common_headers.h
        ./src/code/cgame/cg_camera.h
        ./src/code/Ragl/ragl_common.h
        ./src/code/Ragl/graph_vs.h
        ./src/code/Ragl/graph_region.h
        ./src/code/Ragl/graph_triangulate.h
        ./src/code/Ragl/kdtree_vs.h
        ./src/code/ui/menudef.h
        ./src/code/ui/ui_local.h
        ./src/code/ui/ui_shared.h
        ./src/code/ui/ui_public.h
        ./src/code/ui/gameinfo.h
        ./src/code/rd-common/mdx_format.h
        ./src/code/rd-common/tr_public.h
        ./src/code/rd-common/tr_font.h
        ./src/code/rd-common/tr_types.h
        ./src/code/rd-common/tr_common.h
        ./src/code/Ratl/ratl_common.h
        ./src/code/Ratl/array_vs.h
        ./src/code/Ratl/string_vs.h
        ./src/code/Ratl/list_vs.h
        ./src/code/Ratl/map_vs.h
        ./src/code/Ratl/vector_vs.h
        ./src/code/Ratl/pool_vs.h
        ./src/code/Ratl/hash_pool_vs.h
        ./src/code/Ratl/heap_vs.h
        ./src/code/Ratl/bits_vs.h
        ./src/code/Ratl/handle_pool_vs.h
        ./src/code/Ratl/queue_vs.h
        ./src/code/Ratl/scheduler_vs.h
        ./src/code/Ratl/grid_vs.h
        ./src/code/Ratl/stack_vs.h
        ./src/code/mp3code/l3.h
        ./src/code/mp3code/mp3struct.h
        ./src/code/mp3code/port.h
        ./src/code/mp3code/config.h
        ./src/code/mp3code/small_header.h
        ./src/code/mp3code/tableawd.h
        ./src/code/mp3code/mhead.h
        ./src/code/mp3code/htable.h
        ./src/code/mp3code/jdw.h
        ./src/code/client/vmachine.h
        ./src/code/client/client.h
        ./src/code/client/snd_public.h
        ./src/code/client/snd_local.h
        ./src/code/client/OpenAL/al.h
        ./src/code/client/OpenAL/altypes.h
        ./src/code/client/OpenAL/alu.h
        ./src/code/client/OpenAL/alut.h
        ./src/code/client/OpenAL/alctypes.h
        ./src/code/client/OpenAL/alc.h
        ./src/code/client/keycodes.h
        ./src/code/client/cl_mp3.h
        ./src/code/client/eax/eax.h
        ./src/code/client/eax/EaxMan.h
        ./src/code/client/snd_ambient.h
        ./src/code/client/snd_music.h
        ./src/code/client/client_ui.h
        ./src/code/client/keys.h
        ./src/code/server/server.h
        ./src/code/server/exe_headers.h
        ./src/code/qcommon/stv_version.h
        ./src/code/qcommon/cm_local.h
        ./src/code/qcommon/ojk_saved_game_helper_fwd.h
        ./src/code/qcommon/qfiles.h
        ./src/code/qcommon/ojk_saved_game.h
        ./src/code/qcommon/cm_polylib.h
        ./src/code/qcommon/timing.h
        ./src/code/qcommon/tags.h
        ./src/code/qcommon/qcommon.h
        ./src/code/qcommon/stringed_ingame.h
        ./src/code/qcommon/tri_coll_test.h
        ./src/code/qcommon/stringed_interface.h
        ./src/code/qcommon/chash.h
        ./src/code/qcommon/cm_public.h
        ./src/code/qcommon/hstring.h
        ./src/code/qcommon/MiniHeap.h
        ./src/code/qcommon/ojk_saved_game_helper.h
        ./src/code/qcommon/cm_patch.h
        ./src/code/qcommon/safe/memory.h
        ./src/code/qcommon/matcomp.h
        ./src/code/qcommon/strippublic.h
        ./src/code/qcommon/ojk_saved_game_class_archivers.h
        ./src/code/qcommon/ojk_i_saved_game.h
        ./src/code/qcommon/sstring.h
        ./src/code/qcommon/ojk_scope_guard.h
        ./src/code/qcommon/q_shared.h
        ./src/code/Rufl/hstring.h
        ./src/code/Rufl/hfile.h
        ./src/code/icarus/StdAfx.h
        ./src/code/icarus/sequencer.h
        ./src/code/icarus/IcarusImplementation.h
        ./src/code/icarus/blockstream.h
        ./src/code/icarus/IcarusInterface.h
        ./src/code/icarus/sequence.h
        ./src/code/icarus/taskmanager.h
        ./src/lib/gsl-lite/include/gsl/gsl-lite.h
        ./src/lib/minizip/include/minizip/unzip.h
        ./src/shared/qcommon/q_color.c
        ./src/shared/qcommon/q_string.c
        ./src/shared/qcommon/q_math.c
        ./src/code/mp3code/csbtl3.c
        ./src/code/mp3code/msis.c
        ./src/code/mp3code/uph.c
        ./src/code/mp3code/wavep.c
        ./src/code/mp3code/cwin.c
        ./src/code/mp3code/cdct.c
        ./src/code/mp3code/l3init.c
        ./src/code/mp3code/cwinm.c
        ./src/code/mp3code/cwinb.c
        ./src/code/mp3code/cupl3.c
        ./src/code/mp3code/hwin.c
        ./src/code/mp3code/csbtb.c
        ./src/code/mp3code/mhead.c
        ./src/code/mp3code/towave.c
        ./src/code/mp3code/upsf.c
        ./src/code/mp3code/csbt.c
        ./src/code/mp3code/mdct.c
        ./src/code/mp3code/cupl1.c
        ./src/code/mp3code/cupini.c
        ./src/code/mp3code/cup.c
        ./src/code/mp3code/l3dq.c
        ./src/lib/minizip/ioapi.c
        ./src/lib/minizip/unzip.c
        ./src/shared/sdl/sdl_sound.cpp
        ./src/shared/sdl/sdl_window.cpp
        ./src/shared/sdl/sdl_input.cpp
        ./src/shared/sys/sys_unix.cpp
        ./src/shared/sys/sys_main.cpp
        ./src/shared/sys/con_log.cpp
        ./src/shared/sys/con_tty.cpp
        ./src/shared/sys/snapvector.cpp
        ./src/shared/sys/sys_event.cpp
        ./src/shared/qcommon/safe/string.cpp
        ./src/shared/qcommon/safe/files.cpp
        ./src/code/game/g_breakable.cpp
        ./src/code/game/AI_Howler.cpp
        ./src/code/game/FighterNPC.cpp
        ./src/code/game/wp_trip_mine.cpp
        ./src/code/game/g_objectives.cpp
        ./src/code/game/g_target.cpp
        ./src/code/game/AI_Utils.cpp
        ./src/code/game/g_client.cpp
        ./src/code/game/AI_AssassinDroid.cpp
        ./src/code/game/g_turret.cpp
        ./src/code/game/wp_stun_baton.cpp
        ./src/code/game/g_spawn.cpp
        ./src/code/game/wp_melee.cpp
        ./src/code/game/AI_Sentry.cpp
        ./src/code/game/g_session.cpp
        ./src/code/game/AI_Grenadier.cpp
        ./src/code/game/wp_bot_laser.cpp
        ./src/code/game/NPC_move.cpp
        ./src/code/game/g_fx.cpp
        ./src/code/game/NPC_sounds.cpp
        ./src/code/game/g_misc.cpp
        ./src/code/game/g_weaponLoad.cpp
        ./src/code/game/g_ref.cpp
        ./src/code/game/AI_Droid.cpp
        ./src/code/game/AI_ImperialProbe.cpp
        ./src/code/game/wp_blaster_pistol.cpp
        ./src/code/game/wp_thermal.cpp
        ./src/code/game/wp_noghri_stick.cpp
        ./src/code/game/g_main.cpp
        ./src/code/game/Q3_Interface.cpp
        ./src/code/game/AI_Tusken.cpp
        ./src/code/game/wp_blaster_rifle.cpp
        ./src/code/game/AI_RocketTrooper.cpp
        ./src/code/game/wp_tusken.cpp
        ./src/code/game/g_combat.cpp
        ./src/code/game/AI_Stormtrooper.cpp
        ./src/code/game/AI_BobaFett.cpp
        ./src/code/game/AI_Wampa.cpp
        ./src/code/game/AI_Animal.cpp
        ./src/code/game/AI_HazardTrooper.cpp
        ./src/code/game/AI_Atst.cpp
        ./src/code/game/NPC_behavior.cpp
        ./src/code/game/wp_rocket_launcher.cpp
        ./src/code/game/g_weapon.cpp
        ./src/code/game/g_rail.cpp
        ./src/code/game/g_misc_model.cpp
        ./src/code/game/wp_emplaced_gun.cpp
        ./src/code/game/AI_Mark1.cpp
        ./src/code/game/g_functions.cpp
        ./src/code/game/AI_Remote.cpp
        ./src/code/game/wp_repeater.cpp
        ./src/code/game/AnimalNPC.cpp
        ./src/code/game/AI_SandCreature.cpp
        ./src/code/game/g_nav.cpp
        ./src/code/game/wp_det_pack.cpp
        ./src/code/game/g_savegame.cpp
        ./src/code/game/AI_Jedi.cpp
        ./src/code/game/g_items.cpp
        ./src/code/game/SpeederNPC.cpp
        ./src/code/game/NPC_utils.cpp
        ./src/code/game/AI_Default.cpp
        ./src/code/game/g_roff.cpp
        ./src/code/game/g_object.cpp
        ./src/code/game/wp_concussion.cpp
        ./src/code/game/NPC_misc.cpp
        ./src/code/game/g_camera.cpp
        ./src/code/game/AI_SaberDroid.cpp
        ./src/code/game/AI_Rancor.cpp
        ./src/code/game/g_trigger.cpp
        ./src/code/game/NPC.cpp
        ./src/code/game/g_emplaced.cpp
        ./src/code/game/wp_saber.cpp
        ./src/code/game/AI_Seeker.cpp
        ./src/code/game/genericparser2.cpp
        ./src/code/game/AI_Sniper.cpp
        ./src/code/game/g_active.cpp
        ./src/code/game/g_navigator.cpp
        ./src/code/game/g_itemLoad.cpp
        ./src/code/game/g_navnew.cpp
        ./src/code/game/wp_demp2.cpp
        ./src/code/game/wp_saberLoad.cpp
        ./src/code/game/g_missile.cpp
        ./src/code/game/g_mem.cpp
        ./src/code/game/bg_misc.cpp
        ./src/code/game/bg_slidemove.cpp
        ./src/code/game/g_usable.cpp
        ./src/code/game/wp_atst.cpp
        ./src/code/game/g_vehicles.cpp
        ./src/code/game/WalkerNPC.cpp
        ./src/code/game/NPC_combat.cpp
        ./src/code/game/g_mover.cpp
        ./src/code/game/g_cmds.cpp
        ./src/code/game/bg_pmove.cpp
        ./src/code/game/AI_Mark2.cpp
        ./src/code/game/G_Timer.cpp
        ./src/code/game/NPC_reactions.cpp
        ./src/code/game/bg_vehicleLoad.cpp
        ./src/code/game/NPC_senses.cpp
        ./src/code/game/NPC_stats.cpp
        ./src/code/game/AI_Interrogator.cpp
        ./src/code/game/bg_pangles.cpp
        ./src/code/game/AI_MineMonster.cpp
        ./src/code/game/wp_flechette.cpp
        ./src/code/game/NPC_goal.cpp
        ./src/code/game/bg_panimate.cpp
        ./src/code/game/AI_Civilian.cpp
        ./src/code/game/wp_bowcaster.cpp
        ./src/code/game/g_utils.cpp
        ./src/code/game/NPC_spawn.cpp
        ./src/code/game/AI_GalakMech.cpp
        ./src/code/game/wp_disruptor.cpp
        ./src/code/game/g_inventory.cpp
        ./src/code/game/g_svcmds.cpp
        ./src/code/rd-vanilla/G2_bolts.cpp
        ./src/code/rd-vanilla/tr_shade.cpp
        ./src/code/rd-vanilla/tr_marks.cpp
        ./src/code/rd-vanilla/G2_surfaces.cpp
        ./src/code/rd-vanilla/tr_image.cpp
        ./src/code/rd-vanilla/tr_model.cpp
        ./src/code/rd-vanilla/tr_WorldEffects.cpp
        ./src/code/rd-vanilla/tr_arb.cpp
        ./src/code/rd-vanilla/tr_shadows.cpp
        ./src/code/rd-vanilla/G2_API.cpp
        ./src/code/rd-vanilla/tr_ghoul2.cpp
        ./src/code/rd-vanilla/tr_mesh.cpp
        ./src/code/rd-vanilla/tr_curve.cpp
        ./src/code/rd-vanilla/G2_misc.cpp
        ./src/code/rd-vanilla/tr_stl.cpp
        ./src/code/rd-vanilla/tr_cmds.cpp
        ./src/code/rd-vanilla/tr_subs.cpp
        ./src/code/rd-vanilla/tr_shade_calc.cpp
        ./src/code/rd-vanilla/tr_scene.cpp
        ./src/code/rd-vanilla/tr_quicksprite.cpp
        ./src/code/rd-vanilla/tr_bsp.cpp
        ./src/code/rd-vanilla/tr_light.cpp
        ./src/code/rd-vanilla/tr_shader.cpp
        ./src/code/rd-vanilla/tr_sky.cpp
        ./src/code/rd-vanilla/tr_skin.cpp
        ./src/code/rd-vanilla/tr_world.cpp
        ./src/code/rd-vanilla/tr_surface.cpp
        ./src/code/rd-vanilla/tr_surfacesprites.cpp
        ./src/code/rd-vanilla/tr_main.cpp
        ./src/code/rd-vanilla/tr_draw.cpp
        ./src/code/rd-vanilla/tr_backend.cpp
        ./src/code/rd-vanilla/tr_init.cpp
        ./src/code/rd-vanilla/G2_bones.cpp
        ./src/code/Ravl/CBounds.cpp
        ./src/code/Ravl/CVec.cpp
        ./src/code/cgame/cg_camera.cpp
        ./src/code/cgame/cg_credits.cpp
        ./src/code/cgame/cg_snapshot.cpp
        ./src/code/cgame/FX_HeavyRepeater.cpp
        ./src/code/cgame/FX_NoghriShot.cpp
        ./src/code/cgame/cg_info.cpp
        ./src/code/cgame/cg_marks.cpp
        ./src/code/cgame/cg_servercmds.cpp
        ./src/code/cgame/FX_ATSTMain.cpp
        ./src/code/cgame/cg_event.cpp
        ./src/code/cgame/FX_BryarPistol.cpp
        ./src/code/cgame/cg_main.cpp
        ./src/code/cgame/FX_Emplaced.cpp
        ./src/code/cgame/FX_RocketLauncher.cpp
        ./src/code/cgame/FxSystem.cpp
        ./src/code/cgame/cg_weapons.cpp
        ./src/code/cgame/FxScheduler.cpp
        ./src/code/cgame/cg_syscalls.cpp
        ./src/code/cgame/cg_consolecmds.cpp
        ./src/code/cgame/cg_drawtools.cpp
        ./src/code/cgame/cg_predict.cpp
        ./src/code/cgame/cg_lights.cpp
        ./src/code/cgame/FX_DEMP2.cpp
        ./src/code/cgame/FX_Disruptor.cpp
        ./src/code/cgame/cg_draw.cpp
        ./src/code/cgame/cg_text.cpp
        ./src/code/cgame/FxPrimitives.cpp
        ./src/code/cgame/cg_scoreboard.cpp
        ./src/code/cgame/cg_players.cpp
        ./src/code/cgame/cg_headers.cpp
        ./src/code/cgame/FX_Flechette.cpp
        ./src/code/cgame/cg_localents.cpp
        ./src/code/cgame/cg_playerstate.cpp
        ./src/code/cgame/cg_view.cpp
        ./src/code/cgame/FX_TuskenShot.cpp
        ./src/code/cgame/FX_Bowcaster.cpp
        ./src/code/cgame/cg_effects.cpp
        ./src/code/cgame/cg_ents.cpp
        ./src/code/cgame/FX_Concussion.cpp
        ./src/code/cgame/FxTemplate.cpp
        ./src/code/cgame/FX_Blaster.cpp
        ./src/code/cgame/FxUtil.cpp
        ./src/code/ui/gameinfo.cpp
        ./src/code/ui/ui_syscalls.cpp
        ./src/code/ui/ui_saber.cpp
        ./src/code/ui/ui_shared.cpp
        ./src/code/ui/ui_connect.cpp
        ./src/code/ui/ui_atoms.cpp
        ./src/code/ui/ui_main.cpp
        ./src/code/rd-common/tr_image_load.cpp
        ./src/code/rd-common/tr_noise.cpp
        ./src/code/rd-common/tr_image_png.cpp
        ./src/code/rd-common/tr_image_jpg.cpp
        ./src/code/rd-common/tr_font.cpp
        ./src/code/rd-common/tr_image_tga.cpp
        ./src/code/Ratl/ratl.cpp
        ./src/code/client/cl_scrn.cpp
        ./src/code/client/cl_console.cpp
        ./src/code/client/cl_keys.cpp
        ./src/code/client/cl_cin.cpp
        ./src/code/client/cl_cgame.cpp
        ./src/code/client/cl_ui.cpp
        ./src/code/client/cl_mp3.cpp
        ./src/code/client/snd_dma.cpp
        ./src/code/client/cl_input.cpp
        ./src/code/client/cl_main.cpp
        ./src/code/client/vmachine.cpp
        ./src/code/client/snd_music.cpp
        ./src/code/client/snd_ambient.cpp
        ./src/code/client/snd_mix.cpp
        ./src/code/client/snd_mem.cpp
        ./src/code/client/cl_parse.cpp
        ./src/code/server/sv_ccmds.cpp
        ./src/code/server/sv_snapshot.cpp
        ./src/code/server/exe_headers.cpp
        ./src/code/server/sv_init.cpp
        ./src/code/server/sv_game.cpp
        ./src/code/server/sv_main.cpp
        ./src/code/server/sv_savegame.cpp
        ./src/code/server/sv_world.cpp
        ./src/code/server/sv_client.cpp
        ./src/code/qcommon/stringed_ingame.cpp
        ./src/code/qcommon/stringed_interface.cpp
        ./src/code/qcommon/q_shared.cpp
        ./src/code/qcommon/cm_patch.cpp
        ./src/code/qcommon/ojk_saved_game.cpp
        ./src/code/qcommon/md4.cpp
        ./src/code/qcommon/tri_coll_test.cpp
        ./src/code/qcommon/cm_polylib.cpp
        ./src/code/qcommon/matcomp.cpp
        ./src/code/qcommon/common.cpp
        ./src/code/qcommon/cvar.cpp
        ./src/code/qcommon/cm_test.cpp
        ./src/code/qcommon/cm_trace.cpp
        ./src/code/qcommon/z_memman_pc.cpp
        ./src/code/qcommon/cm_load.cpp
        ./src/code/qcommon/net_chan.cpp
        ./src/code/qcommon/msg.cpp
        ./src/code/qcommon/cmd.cpp
        ./src/code/qcommon/persistence.cpp
        ./src/code/qcommon/files.cpp
        ./src/code/qcommon/hstring.cpp
        ./src/code/qcommon/strip.cpp
        ./src/code/Rufl/hfile.cpp
        ./src/code/Rufl/hstring.cpp
        ./src/code/icarus/TaskManager.cpp
        ./src/code/icarus/BlockStream.cpp
        ./src/code/icarus/Sequence.cpp
        ./src/code/icarus/IcarusImplementation.cpp
        ./src/code/icarus/Sequencer.cpp
        src/globals.cpp src/code/cgame/animtable.cpp)

include("${VITASDK}/share/vita.cmake" REQUIRED)
set(VITA_APP_NAME "StaticJK: Jedi Academy")
set(VITA_TITLEID  "STATICJK1")
set(VITA_VERSION  "01.00")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I./src -flto -I$(VITASDK)/arm-vita-eabi/include/SDL2 -DNDEBUG -DFINAL_BUILD -D__FLOAT_WORD_ORDER=1 -DVITA -D__linux__ -DGIT_VERSION=\"${GIT_HASH}\" -g -O3 -mfpu=neon -fno-optimize-sibling-calls -fsingle-precision-constant -march=armv7-a -mtune=cortex-a9 -mfloat-abi=hard -ffast-math -fsigned-char -fno-short-enums")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I./src -flto -I$(VITASDK)/arm-vita-eabi/include/SDL2 -DNDEBUG -DFINAL_BUILD -D__GNU_VISIBLE -D__FLOAT_WORD_ORDER=1 -DVITA -D__linux__ -DGIT_VERSION=\"${GIT_HASH}\" -g -O3 -mfpu=neon -fno-optimize-sibling-calls -fsingle-precision-constant -march=armv7-a -mtune=cortex-a9 -mfloat-abi=hard -ffast-math -fpermissive -fsigned-char -fno-short-enums")

if (NOT CMAKE_SYSTEM_NAME MATCHES "Emscripten")
    target_link_libraries(${PROJECT_NAME}.elf
		vitagprof
		SDL2
		jpeg
		png
		z
		vitaGL
		vitashark
		SceShaccCgExt
		taihen_stub
		SceDisplay_stub
		SceGxm_stub
		SceIme_stub
		SceAudioIn_stub
		SceSysmodule_stub
		SceCommonDialog_stub
		SceKernelThreadMgr_stub
		SceCtrl_stub
		SceKernelDmacMgr_stub
		SceTouch_stub
		SceAudio_stub
		ScePower_stub
		SceNet_stub
		SceAppMgr_stub
		SceAppUtil_stub
		SceShaccCg_stub
		SceHid_stub
		SceMotion_stub
		c
		stdc++
		m
		mathneon
	)
else ()
    target_link_libraries(${PROJECT_NAME}.elf "-s TOTAL_MEMORY=52428800")
    target_link_libraries(${PROJECT_NAME}.elf "-s ALLOW_MEMORY_GROWTH=1")
    target_link_libraries(${PROJECT_NAME}.elf "-s DEMANGLE_SUPPORT=1")
    target_link_libraries(${PROJECT_NAME}.elf "-s WASM=1")
    target_link_libraries(${PROJECT_NAME}.elf "-s USE_SDL=2")
    target_link_libraries(${PROJECT_NAME}.elf "-s ASSERTIONS=0")
    target_link_libraries(${PROJECT_NAME}.elf "-s SDL2_IMAGE_FORMATS='[\"bmp\",\"png\",\"jpg\"]'")
    target_link_libraries(${PROJECT_NAME}.elf "-s USE_REGAL=1")
    target_link_libraries(${PROJECT_NAME}.elf "-O0")
    target_link_libraries(${PROJECT_NAME}.elf "-s OFFSCREENCANVAS_SUPPORT=1")
    target_link_libraries(${PROJECT_NAME}.elf "-s USE_LIBJPEG=1")
    target_link_libraries(${PROJECT_NAME}.elf "-s USE_LIBPNG=1")
    target_link_libraries(${PROJECT_NAME}.elf "-s FORCE_FILESYSTEM=1")
    target_link_libraries(${PROJECT_NAME}.elf "-s USE_SDL_IMAGE=2")
    target_link_libraries(${PROJECT_NAME}.elf "-s USE_ZLIB=1")
    target_link_libraries(${PROJECT_NAME}.elf "-s ERROR_ON_UNDEFINED_SYMBOLS=0")
    target_link_libraries(${PROJECT_NAME}.elf "-o ${PROJECT_NAME}.elf.html")
    target_link_libraries(${PROJECT_NAME}.elf "-s DEMANGLE_SUPPORT=1")
    target_link_libraries(${PROJECT_NAME}.elf idbfs.js)
    target_compile_options(${PROJECT_NAME}.elf PUBLIC "-O2")
    target_link_libraries(${PROJECT_NAME}.elf "-O2")
    target_link_libraries(${PROJECT_NAME}.elf "-s TOTAL_MEMORY=52428800")
    target_link_libraries(${PROJECT_NAME}.elf "--pre-js ../scripts/preload.js")
    set_target_properties(${PROJECT_NAME}.elf PROPERTIES LINK_FLAGS "-s USE_SDL=2")
endif ()

vita_create_self(${PROJECT_NAME}.self ${PROJECT_NAME}.elf UNSAFE NOASLR)
	
set(VITA_MKSFOEX_FLAGS "-d ATTRIBUTE2=12")
vita_create_vpk(${PROJECT_NAME}.vpk ${VITA_TITLEID} ${PROJECT_NAME}.self
	VERSION ${VITA_VERSION}
	NAME ${VITA_APP_NAME}
	FILE ${PROJECT_SOURCE_DIR}/sce_sys/icon0.png sce_sys/icon0.png
	FILE ${PROJECT_SOURCE_DIR}/sce_sys/livearea/contents/bg.png sce_sys/livearea/contents/bg.png
	FILE ${PROJECT_SOURCE_DIR}/sce_sys/livearea/contents/startup.png sce_sys/livearea/contents/startup.png
	FILE ${PROJECT_SOURCE_DIR}/sce_sys/livearea/contents/template.xml sce_sys/livearea/contents/template.xml
)
