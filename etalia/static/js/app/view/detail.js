define([
    'jquery',
    'app',
    'text!app/templates/detail.html',
    'app/model/detail'
], function ($, App, template) {

    var defaults = {};

    return App.View.Detail = App.Backbone.View.extend({
        tagName: 'div',
        id: 'detail',

        template: App.Handlebars.compile(template),

        events: {
            "click #detail-close": "close"
        },

        initialize: function (options) {
            App.defaults(options, defaults);

            if (!this.model instanceof App.Model.Detail) {
                throw 'Unexpected model type';
            }
        },

        render: function () {

            $('detail-placeholder').replaceWith(this.$el.html(this.template({})));

            this.model
                .get('view').render().$el
                .appendTo(this.$('.document .wrapper'));

            return this.open();
        },

        open: function() {
            $('body').addClass('detail-opened');

            return this;
        },

        close: function() {
            $('body').removeClass('detail-opened');

            return this;
        }
    });
});
