define(['app/app'], function (App) {

    App.Model.DetailButton = App.Backbone.Model.extend({
        _class: 'App.Model.DetailButton',

        defaults: {
            icon: null,
            title: null,
            caption: null,
            callback: null
        },

        validate: function(attrs, options) {
            // TODO
        }
    });

    App.Model.Detail = App.Backbone.Model.extend({
        defaults: {
            view: null,
            left_button: null,
            center_button: null,
            right_button: null
        },

        initialize: function (options) {
            if (!options.view) {
                throw 'options.view is mandatory.';
            }
        },

        _validateButton: function (button) {
            if (0 == String(button.icon).length) {
                throw 'button.icon is mandatory';
            }
            if (0 == String(button.title).length) {
                throw 'button.title is mandatory';
            }
            if (!(typeof button.callback == 'function')) {
                throw 'button.callback is mandatory';
            }
        },

        _setButton: function (path, options) {
            this._validateButton(options);

            this.set(path, new App.Model.DetailButton(options));
        },

        setLeftButton: function (options) {
            options = App._.extend({
                icon: 'left',
                title: null,
                caption: null,
                callback: null
            }, options);

            this._setButton('left_button', options);

            return this;
        },

        setCenterButton: function (options) {
            options = App._.extend({
                icon: 'close',
                title: null,
                caption: null,
                callback: null
            }, options);

            this._setButton('center_button', options);

            return this;
        },

        setRightButton: function (options) {
            options = App._.extend({
                icon: 'right',
                title: null,
                caption: null,
                callback: null
            }, options);

            this._setButton('right_button', options);

            return this;
        }
    });

    return App.Model.Detail;
});
