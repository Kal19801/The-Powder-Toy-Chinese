#pragma once
// Debug-only self test for the EM field (EMWave2 port) and the settings window
// crash regression. Compiled out entirely unless EMFIELD_DEBUG is defined.
// Run the app with a debug build; the state machine below drives itself from
// GameController::Tick and prints [EMSELFTEST] lines to stdout.

class GameModel;
class GameController;

void EMSelfTestTick(GameModel &gameModel, GameController &gameController);
