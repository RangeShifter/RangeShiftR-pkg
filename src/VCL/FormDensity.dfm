object frmDensity: TfrmDensity
  Left = 0
  Top = 0
  Caption = 'Stage-specific density dependence - Stages'#39' weights'
  ClientHeight = 434
  ClientWidth = 851
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -13
  Font.Name = 'Tahoma'
  Font.Style = [fsBold]
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 16
  object Label4: TLabel
    Left = 8
    Top = 8
    Width = 66
    Height = 16
    Caption = 'Fecundity '
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold, fsItalic]
    ParentFont = False
  end
  object Label1: TLabel
    Left = 288
    Top = 8
    Width = 89
    Height = 16
    Caption = 'Development '
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold, fsItalic]
    ParentFont = False
  end
  object Label2: TLabel
    Left = 568
    Top = 8
    Width = 56
    Height = 16
    Caption = 'Survival '
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold, fsItalic]
    ParentFont = False
  end
  object FecundityDens: TStringGrid
    Left = 8
    Top = 30
    Width = 274
    Height = 175
    BiDiMode = bdLeftToRight
    ColCount = 3
    DefaultColWidth = 50
    DoubleBuffered = False
    Enabled = False
    RowCount = 3
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
    ParentBiDiMode = False
    ParentDoubleBuffered = False
    ParentFont = False
    TabOrder = 0
    Visible = False
  end
  object DevelopDens: TStringGrid
    Left = 288
    Top = 30
    Width = 274
    Height = 175
    BiDiMode = bdLeftToRight
    ColCount = 3
    DefaultColWidth = 50
    DoubleBuffered = False
    Enabled = False
    RowCount = 3
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
    ParentBiDiMode = False
    ParentDoubleBuffered = False
    ParentFont = False
    TabOrder = 1
    Visible = False
  end
  object SurvivalDens: TStringGrid
    Left = 568
    Top = 30
    Width = 274
    Height = 175
    BiDiMode = bdLeftToRight
    ColCount = 3
    DefaultColWidth = 50
    DoubleBuffered = False
    Enabled = False
    RowCount = 3
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    Options = [goFixedVertLine, goFixedHorzLine, goVertLine, goHorzLine, goRangeSelect, goEditing, goTabs, goAlwaysShowEditor]
    ParentBiDiMode = False
    ParentDoubleBuffered = False
    ParentFont = False
    TabOrder = 2
    Visible = False
  end
  object BtnOK: TButton
    Left = 658
    Top = 229
    Width = 89
    Height = 25
    Caption = 'OK'
    Default = True
    TabOrder = 3
    OnClick = BtnOKClick
  end
  object BtnCancel: TButton
    Left = 753
    Top = 229
    Width = 89
    Height = 25
    Cancel = True
    Caption = 'Cancel'
    ModalResult = 2
    TabOrder = 4
    OnClick = BtnCancelClick
  end
end
