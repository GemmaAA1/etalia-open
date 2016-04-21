define(['jquery', 'app', 'app/model/user/user-lib'], function ($, App) {

    App.Model.User = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/user/users',
        //urlRoot: 'http://sf-demo.jessie.dev/app_dev.php/user/users',

        defaults: {
            email: null,
            first_name: null,
            last_name: null,
            photo_url: null
        },

        relations: [
            {
                type: App.Backbone.HasOne,
                key: 'user_lib',
                relatedModel: App.Model.UserLib,
                includeInJSON: 'link',
                autoFetch: false,
                reverseRelation: {
                    key: 'user',
                    type: App.Backbone.HasOne,
                    includeInJson: 'link'
                }
            }
        ],

        schema: {
            email: {type: 'Text', validators: ['required']},
            first_name: {type: 'Text', validators: ['required']},
            last_name: {type: 'Text', validators: ['required']}
            //photo_url: {type: 'Text', validators: ['required']}
        }
    });

    var currentUser = null;

    /**
     * Returns the current (authenticated) user.
     *
     * @return UserModel
     *
     * @TODO fetch current (authenticated) user id
     */
    App.Model.User.getCurrent = function() {
        if (!currentUser) {
            var id = $('body').data('user-id');
            if (!id) {
                // TODO throw 'Failed to determine user id';
                id = 21;
            }
            currentUser = App.Model.User.find(id);
            if (!currentUser) {
                currentUser = new App.Model.User({id: id});
                currentUser.fetch();
            }
        }
        return currentUser;
    };

    return App.Model.User;
});
